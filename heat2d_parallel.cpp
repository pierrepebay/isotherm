#include <mpi.h>

#ifdef VTK_AVAILABLE
#include <vtkNamedColors.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageMapToColors.h>
#include <vtkLookupTable.h>
#include <vtkPNGWriter.h>
#endif

#include <iostream>
#include <vector>
#include <array>
#include <numeric>
#include <algorithm>
#include <fstream>

// Utility function to convert 2D coordinates to a single index in a linear array
uint64_t cartesianToIndex(uint64_t i, uint64_t j, uint64_t N_local) {
  return j * N_local + i;
}

void exportDataToText(const std::vector<double> &T, const uint64_t iteration, uint64_t N_global) {
  std::ofstream outFile("output/temperature_data_" + std::to_string(iteration) + ".txt");
  for (uint64_t j = 0; j < N_global; ++j) {
    for (uint64_t i = 0; i < N_global; ++i) {
      outFile << T[j * N_global + i];
      if (i != N_global - 1) outFile << ", ";
    }
    outFile << "\n";
  }
  outFile.close();
}

#ifdef VTK_AVAILABLE
void renderData(const std::vector<double> &T, const uint64_t iteration, uint64_t N_global) {
  vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
  imageData->SetDimensions(N_global, N_global, 1);
  imageData->AllocateScalars(VTK_DOUBLE, 1);

  // Fill image data with temperature values
  for (int j = 0; j < N_global; ++j) {
    for (int i = 0; i < N_global; ++i) {
      double* pixel = static_cast<double*>(imageData->GetScalarPointer(i, j, 0));
      pixel[0] = T[j * N_global + i];
    }
  }

  // Create a lookup table to map scalar values to colors
  vtkSmartPointer<vtkLookupTable> lookupTable = vtkSmartPointer<vtkLookupTable>::New();
  lookupTable->SetNumberOfTableValues(256); // 256 colors
  lookupTable->SetRange(0.0, 100.0); // Temperature range
  lookupTable->Build();

  // Fill the lookup table with colors
  for (int i = 0; i < 256; ++i) {
    double r = i / 255.0; // gradient from cool (blue) to warm (red)
    double g = (255 - i) / 255.0;
    double b = 1.0 - r;
    lookupTable->SetTableValue(i, r, g, b, 1.0); // RGBA, A=1.0 means fully opaque
  }

  // Map the scalar values in the image to colors using the lookup table
  vtkSmartPointer<vtkImageMapToColors> mapColors = vtkSmartPointer<vtkImageMapToColors>::New();
  mapColors->SetLookupTable(lookupTable);
  mapColors->SetInputData(imageData);
  mapColors->Update();

  // Write the color-mapped image data to a PNG file
  vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
  writer->SetFileName(("output/temperature_" + std::to_string(iteration) + ".png").c_str());
  writer->SetInputData(mapColors->GetOutput());
  writer->Write();
}
#endif

#ifdef VTK_AVAILABLE
// Function to export mesh data to a VTK file, used by process 0
void exportMesh(const std::vector<double> &T, const uint64_t iteration, uint64_t N_global) {
  std::cout << "Exporting mesh for iteration " << iteration << " with N_global = " << N_global << std::endl;

  vtkNew<vtkNamedColors> colors;
  vtkNew<vtkImageData> mesh;

  mesh->SetDimensions(N_global + 1, N_global + 1, 1);
  mesh->SetSpacing(1.0, 1.0, 1.0);
  mesh->SetOrigin(0.0, 0.0, 0.0);
  mesh->AllocateScalars(VTK_DOUBLE, 1);

  vtkNew<vtkDoubleArray> scalar_double_array;
  scalar_double_array->SetNumberOfComponents(1);
  scalar_double_array->SetNumberOfTuples(N_global * N_global);
  scalar_double_array->SetName("Temperature");

  for (uint64_t j = 0; j < N_global; ++j) {
    for (uint64_t i = 0; i < N_global; ++i) {
      scalar_double_array->SetTuple1(cartesianToIndex(i, j, N_global), T[cartesianToIndex(i, j, N_global)]);
      // print values of T if iteration is 20
      if (iteration == 20) {
        std::cout << "T[" << i << ", " << j << "] = " << T[cartesianToIndex(i, j, N_global)] << std::endl;
      }
    }
  }

  mesh->GetCellData()->SetScalars(scalar_double_array);

  vtkNew<vtkXMLImageDataWriter> writer;

  // Create output directory if it does not exist
  std::string command = "mkdir -p output";
  system(command.c_str());

  std::string mesh_filename = "output/output" + std::to_string(iteration) + ".vti";
  writer->SetFileName(mesh_filename.c_str());
  writer->SetInputData(mesh);
  writer->Write();
}
#endif

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <grid size>\n";
    return 1;
  }

  const int N = std::stoi(argv[1]);

  double start_time, end_time;

  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Barrier(MPI_COMM_WORLD);
  start_time = MPI_Wtime();

  const int N_local = N / size; // Assuming N is divisible by size
  std::vector<double> T_local(N_local * N, 100.0); // Initialize local temperature array
  std::vector<double> Tnew_local(N_local * N);

  double cx = 0.1, cy = 0.1;
  uint64_t max_iter = 100;

  for (uint64_t iter = 0; iter < max_iter; ++iter) {
    // Prepare data for sending and receiving boundary rows
    std::vector<double> send_up(N, 0), recv_up(N, 0);
    std::vector<double> send_down(N, 0), recv_down(N, 0);

    if (rank > 0) { // Send the first row up and receive from above
      std::copy(T_local.begin(), T_local.begin() + N, send_up.begin());
      MPI_Sendrecv(send_up.data(), N, MPI_DOUBLE, rank - 1, 0,
                   recv_up.data(), N, MPI_DOUBLE, rank - 1, 0,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      std::copy(recv_up.begin(), recv_up.end(), T_local.begin()); // Copy received data to local array
    }

    if (rank < size - 1) { // Send the last row down and receive from below
      std::copy(T_local.end() - N, T_local.end(), send_down.begin());
      MPI_Sendrecv(send_down.data(), N, MPI_DOUBLE, rank + 1, 0,
                   recv_down.data(), N, MPI_DOUBLE, rank + 1, 0,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      std::copy(recv_down.begin(), recv_down.end(), T_local.end() - N); // Copy received data to local array
    }

    // Compute new temperatures for inner cells of local subdomain
    // Ensure not to update the first and last columns, treating them as boundaries
    for (uint64_t j = 1; j < N_local - 1; ++j) { // Skip the first and last rows as they are boundaries
      for (uint64_t i = 1; i < N - 1; ++i) { // Only update the inner cells, not the first and last columns
        Tnew_local[cartesianToIndex(i, j, N)] = T_local[cartesianToIndex(i, j, N)] +
                                                cx * (T_local[cartesianToIndex(i + 1, j, N)] - 2 * T_local[cartesianToIndex(i, j, N)] + T_local[cartesianToIndex(i - 1, j, N)]) +
                                                cy * (T_local[cartesianToIndex(i, j + 1, N)] - 2 * T_local[cartesianToIndex(i, j, N)] + T_local[cartesianToIndex(i, j - 1, N)]);
      }
    }


    // Ensure boundary conditions are maintained
    for (uint64_t j = 0; j < N_local; ++j) {
      Tnew_local[cartesianToIndex(0, j, N)] = 0; // Left boundary of the sub-domain
      Tnew_local[cartesianToIndex(N - 1, j, N)] = 0; // Right boundary of the sub-domain
    }
    if (rank == 0) {
      std::fill(Tnew_local.begin(), Tnew_local.begin() + N, 0); // Top boundary
    }
    if (rank == size - 1) {
      std::fill(Tnew_local.end() - N, Tnew_local.end(), 0); // Bottom boundary
    }

    T_local.swap(Tnew_local); // Update local temperatures

    // Gather all subdomains at process 0
    std::vector<double> T_global;
    if (rank == 0) {
      T_global.resize(N * N);
    }
    MPI_Gather(T_local.data(), N_local * N, MPI_DOUBLE, T_global.data(), N_local * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//    // Process 0 exports the mesh
//    if (rank == 0)
//      exportMesh(T_global, iter, N);
//    }

    // process 0 prints temperature min, max, and average
    if (rank == 0) {
      double T_min = *std::min_element(T_global.begin(), T_global.end());
      double T_max = *std::max_element(T_global.begin(), T_global.end());
      double T_avg = std::accumulate(T_global.begin(), T_global.end(), 0.0) / (N * N);
      std::cout << "Iteration " << iter << ": T_min = " << T_min << ", T_max = " << T_max << ", T_avg = " << T_avg << std::endl;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  end_time = MPI_Wtime();

  MPI_Finalize();

  double execution_time = end_time - start_time;
  if (rank == 0) {
    std::cout << "Execution time: " << execution_time << " seconds\n";
  }

  std::ofstream outFile;
  outFile.open("execution_times.csv", std::ios_base::app); // Append mode
  outFile << N << "," << size << "," << execution_time << "\n";

  return 0;
}
