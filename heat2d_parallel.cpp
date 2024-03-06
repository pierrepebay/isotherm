#include <vtkNamedColors.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>

#include <mpi.h>

#include <iostream>
#include <array>
#include <vector>

const int N = 20;

uint64_t cartesianToIndex(uint64_t i, uint64_t j) {
  return j * N + i;
}

void exportMesh(const std::vector<double> & T, const uint64_t  iteration) {
  vtkNew<vtkNamedColors> colors;
  vtkNew<vtkImageData> mesh;

  mesh->SetDimensions(N + 1, N + 1, 1);
  mesh->SetSpacing(1.0, 1.0, 1.0);
  mesh->SetOrigin(0.0, 0.0, 0.0);
  mesh->AllocateScalars(VTK_DOUBLE, 1);

  vtkNew<vtkDoubleArray> scalar_double_array;
  scalar_double_array->SetNumberOfComponents(1);
  scalar_double_array->SetNumberOfTuples(N * N);
  scalar_double_array->SetName("Temperature");

  for (vtkIdType i = 0; i < N * N; ++i) {
    scalar_double_array->SetValue(i, T[i]);
  }

  mesh->GetCellData()->SetScalars(scalar_double_array);

  vtkNew<vtkXMLImageDataWriter> writer;
  std::string mesh_filneame ="output" + std::to_string(iteration) + ".vti";
  writer->SetFileName(mesh_filneame.c_str());
  writer->SetInputData(mesh);
  writer->Write();
}

void printMesh(const std::vector<double> & T) {
  for (uint64_t i = 0; i < N; ++i) {
    for (uint64_t j = 0; j < N; ++j) {
      std::cout << T[cartesianToIndex(i, j)] << " ";
    }
    std::cout << std::endl;
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  std::cout << "Hello, World!" << std::endl;

  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  // Print off a hello world message
  printf("Hello world from processor %s, rank %d out of %d processors\n",
         processor_name, rank, size);

  int slice_size = N / size;
  std::vector<double> T(slice_size*N);

  double cx = 0.1;
  double cy = 0.1;

  for (double & u : T)
    u = 0.0;

  uint64_t max_iter = 100;

  for (uint64_t iter = 0; iter < max_iter; ++iter) {
    std::vector<double> Tnew(slice_size*N);
    for (uint64_t i = 0; i < slice_size; ++i) {
      for (uint64_t j = 0; j < N; ++j) {
        if (i == 0 || i == slice_size - 1 || j == 0 || j == N - 1) {
          Tnew[cartesianToIndex(i, j)] = 1000.0;
        } else {
          Tnew[cartesianToIndex(i, j)] = T[cartesianToIndex(i, j)] +
                                         cx * (T[cartesianToIndex(i + 1, j)] - 2 * T[cartesianToIndex(i, j)] + T[cartesianToIndex(i - 1, j)]) +
                                         cy * (T[cartesianToIndex(i, j + 1)] - 2 * T[cartesianToIndex(i, j)] + T[cartesianToIndex(i, j - 1)]);
        }
      }
    }

    // Exchange boundary data with neighbors
    if (rank > 0) {
      MPI_Send(&Tnew[0], N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
      MPI_Recv(&Tnew[N*(slice_size-1)], N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (rank < size - 1) {
      MPI_Send(&Tnew[N*(slice_size-1)], N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
      MPI_Recv(&Tnew[0], N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    T = Tnew;

    // Gather results to root process for writing the output
    std::vector<double> T_global(N*N);
    MPI_Gather(&T[0], slice_size*N, MPI_DOUBLE, &T_global[0], slice_size*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
      exportMesh(T_global, iter);
    }
  }

  MPI_Finalize();
  return 0;
}