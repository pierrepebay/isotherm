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
#endif

#include <iostream>
#include <array>

const int N = 20;

uint64_t cartesianToIndex(uint64_t i, uint64_t j) {
  return j * N + i;
}

#ifdef VTK_AVAILABLE
void exportMesh(const std::array<double, N*N> & T, const uint64_t  iteration) {
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

  // create output directory if it does not exist
  std::string command = "mkdir -p output";
  system(command.c_str());

  std::string mesh_filename ="output/output" + std::to_string(iteration) + ".vti";
  writer->SetFileName(mesh_filename.c_str());
  writer->SetInputData(mesh);
  writer->Write();
}
#endif

void printMesh(const std::array<double, N*N> & T) {
  for (uint64_t i = 0; i < N; ++i) {
    for (uint64_t j = 0; j < N; ++j) {
      std::cout << T[cartesianToIndex(i, j)] << " ";
    }
    std::cout << std::endl;
  }
}

int main() {
  std::array<double, N*N> T{};

  double cx = 0.1;
  double cy = 0.1;

  for (double & u : T)
    u = 1000.0;

  uint64_t max_iter = 100;

  for (uint64_t iter = 0; iter < max_iter; ++iter) {
    std::array<double, N*N> Tnew{};
    for (uint64_t j = 0; j < N; ++j) {
      for (uint64_t i = 0; i < N; ++i) {
        if (i == 0 || i == N - 1 || j == 0 || j == N - 1) {
          Tnew[cartesianToIndex(i, j)] = 0.0;
        } else {
          Tnew[cartesianToIndex(i, j)] = T[cartesianToIndex(i, j)] +
            cx * (T[cartesianToIndex(i + 1, j)] - 2 * T[cartesianToIndex(i, j)] + T[cartesianToIndex(i - 1, j)]) +
            cy * (T[cartesianToIndex(i, j + 1)] - 2 * T[cartesianToIndex(i, j)] + T[cartesianToIndex(i, j - 1)]);
        }
      }
    }
    T = Tnew;
//    exportMesh(T, iter);
  }
}
