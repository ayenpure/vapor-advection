#include <chrono>
#include <iostream>
#include <vector>
#include <random>

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

#include <vapor/DataMgr.h>
#include <vapor/FileUtils.h>
#include <vapor/VaporField.h>

#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>

#include <vtkm/Types.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/Invoker.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/worklet/LagrangianStructures.h>
#include <vtkm/worklet/lcs/GridMetaData.h>
#include <vtkm/worklet/WorkletMapTopology.h>

#include "Advection.h"
#include "FTLEHelper.h"

using namespace VAPoR;

int GetFilesFromDirectory(const std::string& datapath,
                           std::vector<std::string>& datafiles)
{
  if ( !boost::filesystem::exists( datapath ) ) return -1;
  boost::filesystem::directory_iterator end_itr; // default construction yields past-the-end
  for ( boost::filesystem::directory_iterator itr( datapath );
        itr != end_itr;
        ++itr )
  {
    datafiles.push_back(itr->path().c_str());
  }
  return 0;
}

namespace detail
{

void VTKmToVaporParticle(const std::vector<vtkm::Vec3f>& vtkm,
                         std::vector<flow::Particle>& vapor)
{
  for(long i = 0; i < vtkm.size(); i++ )
  {
    flow::Particle particle;
    vtkm::Vec3f vec = vtkm.at(i);
    std::cout << "[Copying " << i << "]" << vec << std::endl;
    particle.location.x = vec[0];
    particle.location.y = vec[1];
    particle.location.z = vec[2];
    particle.time       = 0;
    vapor.push_back(particle);
  }
}

class VaporToVTKm : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_EXEC_CONT VaporToVTKm() = default;

  using ControlSignature = void(FieldIn, FieldOut);

  VTKM_EXEC void operator()(const flow::Particle& vaporLoc,
                            vtkm::Vec3f& vtkmLoc) const
  {
    auto x = vaporLoc.location.x;
    auto y = vaporLoc.location.y;
    auto z = vaporLoc.location.z;
    vtkmLoc = vtkm::Vec3f{x, y, z};
  }
};

class CompareWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_EXEC_CONT CompareWorklet() = default;

  using ControlSignature = void(FieldIn, FieldIn);
  using ExecutionSignature = void(WorkIndex, _1, _2);

  VTKM_EXEC void operator()(const vtkm::Id index,
                            const double& vapor,
                            const double& visit) const
  {
    std::cout << "[" << index << "] "<<" VaPOR : " << vapor << " | VTK-m : " << visit << std::endl;
  }
};

/*class GetFlowMap : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
  VTKM_EXEC_CONT GetFlowMap() = default;

  using ControlSignature = void(CellSetIn, WholeArrayIn, FieldOut, FieldOut);
  using ExecutionSignature = void(CellShape, PointCount, PointIndices, _2, _3, _4);

  template <typename CellShapeTag,
            typename PointCountType,
            typename PointIndicesType,
            typename CoordinatesPortalType,
            typename PointType>
  VTKM_EXEC void operator()(CellShapeTag shape,
                            PointCountType pointCount,
                            const PointIndicesType &pointIndices,
                            const CoordinatesPortalType &coordsPortal,
                            PointType& start,
                            PointType& end) const
  {
    start = coordsPortal.Get(pointIndices[0]);
    end   = coordsPortal.Get(pointIndices[pointCount - 1]);
  }
};*/

} //namespace detail

void CopyToVTKmArrays(std::vector<flow::Particle>& _vaporLocation,
                      vtkm::cont::ArrayHandle<vtkm::Vec3f>& vtkmLocation)
{
  vtkm::cont::Invoker invoker;
  vtkm::cont::ArrayHandle<flow::Particle> vaporLocation;
  vaporLocation = vtkm::cont::make_ArrayHandle(_vaporLocation);
  invoker(detail::VaporToVTKm{}, vaporLocation, vtkmLocation);
}

int main (int argc, char** argv)
{
  // Input params
  // 1. File
  // 2. Field variable
  // 3. Number of steps
  // 4. Step length (deltaT) : step length is calculated dynamically
  //    This is suggestive?
  // 5. Seeding parameters
  namespace options = boost::program_options;
  options::options_description desc("Options");
  desc.add_options()("data", options::value<std::string>()->required(), "Mesh data")
                    ("flow", options::value<std::string>()->required(), "Trajectory data");
  options::variables_map vm;
  options::store(options::parse_command_line(argc, argv, desc), vm); // can throw
  options::notify(vm);
  if (!(vm.count("data")
      && vm.count("flow")))
  {
    std::cout << "Advection Benchmark" << std::endl << desc << std::endl;
  }

  std::string datapath = vm["data"].as<std::string>();
  std::string flowpath = vm["flow"].as<std::string>();
  // Field is not really used
  std::cout << "Advection w/ : "
            << "\nData : " << datapath << std::endl;

  std::vector<vtkm::Vec3f> startPoints, endPoints;

  vtkm::cont::DataSet data;
  {
    vtkm::io::VTKDataSetReader reader(datapath);
    data = reader.ReadDataSet();
  }
  {
    vtkSmartPointer<vtkPolyDataReader> inputReader
      = vtkSmartPointer<vtkPolyDataReader>::New();
    inputReader->SetFileName(flowpath.c_str());
    inputReader->Update();

    vtkSmartPointer<vtkPolyData> flowDataset = inputReader->GetOutput();
    vtkSmartPointer<vtkPoints> flowPoints = flowDataset->GetPoints();
    vtkSmartPointer<vtkCellArray> flowCells = flowDataset->GetLines();

    // Get basic data needed for performing the merge.
    vtkIdType  numCells;
    numCells = flowDataset->GetNumberOfCells();
    if(numCells == 0)
    {
      std::cout << "Number of cells is 0, something went wrong" << std::endl;
      exit(EXIT_FAILURE);
    }
    vtkIdType numPoints;
    vtkIdType* points;
    while(flowCells->GetNextCell(numPoints, points))
    {
      vtkIdType startIndex = points[0];
      vtkIdType endIndex = points[numPoints - 1];
      double* startPoint = flowPoints->GetPoint(startIndex);
      vtkm::Vec3f _1 = vtkm::Vec3f{startPoint[0], startPoint[1], startPoint[2]};
      double* endPoint = flowPoints->GetPoint(endIndex);
      vtkm::Vec3f _2 = vtkm::Vec3f{endPoint[0], endPoint[1], endPoint[2]};
      std::cout << "{" << _1 << " : " << _2 << "} : " << (_1 == _2) << std::endl;
      startPoints.push_back(_1);
      endPoints.push_back(_2);
    }
  }

  std::cout << "Number of point : " << startPoints.size() << " : " << endPoints.size() << std::endl;

  vtkm::cont::CoordinateSystem coords = data.GetCoordinateSystem();
  vtkm::cont::DynamicCellSet cellset  = data.GetCellSet();

  std::vector<size_t> dims;
  using Structured2DType = vtkm::cont::CellSetStructured<2>;
  using Structured3DType = vtkm::cont::CellSetStructured<3>;
  if (cellset.IsType<Structured2DType>())
  {
    vtkm::Id2 _dims =
        cellset.Cast<Structured2DType>().GetSchedulingRange(vtkm::TopologyElementTagPoint());
    dims.push_back(_dims[0]);
    dims.push_back(_dims[1]);
    dims.push_back(size_t(1));
  }
  else
  {
    vtkm::Id3 _dims =
        cellset.Cast<Structured3DType>().GetSchedulingRange(vtkm::TopologyElementTagPoint());
    dims.push_back(_dims[0]);
    dims.push_back(_dims[1]);
    dims.push_back(_dims[2]);
  }
  std::vector<double> bounds;
  vtkm::Bounds _bounds = coords.GetBounds();
  bounds.push_back(_bounds.X.Min);
  bounds.push_back(_bounds.X.Max);
  bounds.push_back(_bounds.Y.Min);
  bounds.push_back(_bounds.Y.Max);
  bounds.push_back(_bounds.Z.Min);
  bounds.push_back(_bounds.Z.Max);

  std::cout << "Calculating FTLE using VAPoR" << std::endl;
  std::cout << dims.size() << " | " << bounds.size() << std::endl;

  detail::GridMetaData metaData(dims, bounds);
  std::vector<double> _vaporFTLE;
  std::vector<flow::Particle> vaporStart, vaporEnd;
  detail::VTKmToVaporParticle(startPoints, vaporStart);
  detail::VTKmToVaporParticle(endPoints, vaporEnd);
  CalculateFTLE(vaporStart, vaporEnd, metaData, 10, _vaporFTLE);
  vtkm::cont::ArrayHandle<double> vaporFTLE;
  vaporFTLE = vtkm::cont::make_ArrayHandle(_vaporFTLE);

  std::cout << "Calculating FTLE using VTK-m" << std::endl;

  vtkm::cont::ArrayHandle<double> vtkmFTLE;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> vtkmStart, vtkmEnd;
  vtkmStart = vtkm::cont::make_ArrayHandle(startPoints);
  vtkmEnd = vtkm::cont::make_ArrayHandle(endPoints);
  vtkm::cont::Invoker invoker;
  vtkm::worklet::LagrangianStructures<2> lcsWorklet(10, data.GetCellSet());
  invoker(lcsWorklet, vtkmStart, vtkmEnd, vtkmFTLE);

  std::cout << "Comparing values" << std::endl;
  invoker(detail::CompareWorklet{}, vaporFTLE, vtkmFTLE);

  return 0;
}
