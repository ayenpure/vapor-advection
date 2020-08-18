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
#include <vtkm/Particle.h>

#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/Invoker.h>

#include <vtkm/io/VTKDataSetReader.h>

#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkm/worklet/particleadvection/Field.h>
#include <vtkm/worklet/particleadvection/GridEvaluators.h>
#include <vtkm/worklet/particleadvection/Integrators.h>
#include <vtkm/worklet/particleadvection/Particles.h>

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

class ExtractParticlePosition : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn particle, FieldOut position);
  using ExecutionSignature = void(_1, _2);
  using InputDomain = _1;

  VTKM_EXEC void operator()(const vtkm::Massless& particle, vtkm::Vec3f& pt) const
  {
    pt = particle.Pos;
  }
};

class MakeParticles : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn seed, FieldOut particle);
  using ExecutionSignature = void(WorkIndex, _1, _2);
  using InputDomain = _1;

  VTKM_EXEC void operator()(const vtkm::Id index,
                            const vtkm::Vec3f& seed,
                            vtkm::Massless& particle) const
  {
    particle.ID = index;
    particle.Pos = seed;
  }
};

void VTKmToVaporParticle(const vtkm::cont::ArrayHandle<vtkm::Vec3f>& vtkm,
                         std::vector<flow::Particle>& vapor)
{
  auto portal = vtkm.GetPortalConstControl();

  for(long i = 0; i < vtkm.GetNumberOfValues(); i++ )
  {
    flow::Particle particle;
    vtkm::Vec3f vec = portal.Get(i);
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

  using ControlSignature = void(FieldIn, FieldIn, FieldIn);
  using ExecutionSignature = void(WorkIndex, _1, _2, _3);

  VTKM_EXEC void operator()(const vtkm::Id index,
                            const double& vapor,
                            const double& vtkm,
                            const double& visit) const
  {
    std::cout << "[" << index << "]"
              << " VaPOR : " << vapor
              << " | VTK-m : " << vtkm
              << " | VisIt : " << visit
              << std::endl;
  }
};

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
  desc.add_options()("data", options::value<std::string>()->required(), "Path to dataset")
                    ("field", options::value<std::string>()->required(), "Name of vector field")
                    ("steps", options::value<long>()->required(), "Number of Steps")
                    ("length", options::value<float>()->required(), "Length of a single step");
  options::variables_map vm;
  options::store(options::parse_command_line(argc, argv, desc), vm); // can throw
  options::notify(vm);
  if (!(vm.count("data")
      && vm.count("field")
      && vm.count("steps")
      && vm.count("length")))
  {
    std::cout << "Advection Benchmark" << std::endl << desc << std::endl;
  }

  std::string datapath = vm["data"].as<std::string>();
  // Field is not really used
  std::string fieldname = vm["field"].as<std::string>();
  long steps = vm["steps"].as<long>();
  float length = vm["length"].as<float>();

  vtkm::cont::DataSet data;
  {
    vtkm::io::VTKDataSetReader reader(datapath);
    data = reader.ReadDataSet();
  }
  vtkm::cont::ArrayHandle<double> vtkmFTLE;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> vtkmStart, vtkmEnd;
  vtkm::cont::Invoker invoker;
  {
    using FieldHandle = vtkm::cont::ArrayHandle<vtkm::Vec3f>;
    using FieldType = vtkm::worklet::particleadvection::VelocityField<FieldHandle>;
    using GridEvaluator = vtkm::worklet::particleadvection::GridEvaluator<FieldType>;
    using Integrator = vtkm::worklet::particleadvection::RK4Integrator<GridEvaluator>;

    vtkm::cont::ArrayCopy(data.GetCoordinateSystem().GetData(), vtkmStart);

    FieldHandle field;
    data.GetField(fieldname).GetData().CopyTo(field);
    FieldType velocities(field);
    GridEvaluator evaluator(data.GetCoordinateSystem(), data.GetCellSet(), velocities);
    Integrator integrator(evaluator, length);
    vtkm::worklet::ParticleAdvection particles;
    vtkm::worklet::ParticleAdvectionResult<vtkm::Massless> advectionResult;
    vtkm::cont::ArrayHandle<vtkm::Massless> advectionPoints;
    invoker(detail::MakeParticles{}, vtkmStart, advectionPoints);
    advectionResult = particles.Run(integrator, advectionPoints, steps);
    invoker(detail::ExtractParticlePosition{}, advectionResult.Particles, vtkmEnd);
    // Advect and get start and end points for the trajectories.
  }
  std::cout << "Calculating FTLE using VTK-m" << std::endl;
  vtkm::worklet::LagrangianStructures<2> lcsWorklet(10, data.GetCellSet());
  invoker(lcsWorklet, vtkmStart, vtkmEnd, vtkmFTLE);

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
  detail::VTKmToVaporParticle(vtkmStart, vaporStart);
  detail::VTKmToVaporParticle(vtkmEnd, vaporEnd);
  CalculateFTLE(vaporStart, vaporEnd, metaData, 10, _vaporFTLE);
  vtkm::cont::ArrayHandle<double> vaporFTLE;
  vaporFTLE = vtkm::cont::make_ArrayHandle(_vaporFTLE);

  std::cout << "Comparing values" << std::endl;
  vtkm::cont::ArrayHandle<double> visitFTLE;
  data.GetField("ftle").GetData().CopyTo(visitFTLE);
  invoker(detail::CompareWorklet{}, vaporFTLE, vtkmFTLE, visitFTLE);

  return 0;
}
