#include <chrono>
#include <iostream>
#include <vector>
#include <random>
#include <omp.h>

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

#include "Advection.h"

#include <vapor/DataMgr.h>
#include <vapor/FileUtils.h>
#include <vapor/VaporField.h>

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

void GenerateSeeds(std::vector<flow::Particle>& seeds,
                   const std::vector<double>& xrake,
                   const std::vector<double>& yrake,
                   const std::vector<double>& zrake,
                   const long numOfSeeds,
                   const double initTime)
{
  std::cout << "[begin] generate seeds" << std::endl;
  const unsigned int randSeed = 32;
  std::mt19937 gen(randSeed); //Standard mersenne_twister_engine
  std::uniform_real_distribution<float> distX( xrake.at(0), xrake.at(1));
  std::uniform_real_distribution<float> distY( yrake.at(0), yrake.at(1));
  std::uniform_real_distribution<float> distZ( zrake.at(0), zrake.at(1));
  seeds.resize(1);
  for( long i = 0; i < 1; i++ )
  {
            seeds[i].location.x = 1.448671;
            seeds[i].location.y = -8.197534;
            seeds[i].location.z = 0.344849;
    //seeds[i].location.x = distX(gen);
    //seeds[i].location.y = distY(gen);
    //seeds[i].location.z = distZ(gen);
    seeds[i].time       = initTime;
  }
  for(auto& seed : seeds)
  {
    std::cout << "{" << seed.location.x << "," << seed.location.y << ", " << seed.location.z << "} : " << seed.time
              << std::endl;
  }
  std::cout << "[end] generate seeds" << std::endl;
}

int PrintStreams(external::Advection& advector)
{
  size_t numStreams = advector.GetNumberOfStreams();
  for(size_t index = 0; index < numStreams; index++)
  {
    auto& stream = advector.GetStreamAt(index);
    {
      auto& start = stream.front();
      auto& end = stream.back();
      std::cout << "[" << index << "] Size : " << stream.size() << " "
                << "{(" << start.location.x << ", " << start.location.y << ", " << start.location.z << ") : "
                << "(" << end.location.x << ", " << end.location.y << ", " << end.location.z << ")}" << std::endl;
    }
  }
  return 0;
}

void multitoken_double(const std::string& tokens,
                       std::vector<double>& doubles)
{
  const std::string delimiter = " ";
  size_t pos = 0;
  pos = tokens.find(delimiter);
  doubles.push_back(stod(tokens.substr(0, pos)));
  doubles.push_back(stod(tokens.substr(pos, tokens.size())));
  std::cout << doubles.size() << std::endl;
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
  desc.add_options()("datapath", options::value<std::string>()->required(), "Path to dataset")
                    ("fieldx", options::value<std::string>()->required(), "Name of vector field")
                    ("fieldy", options::value<std::string>()->required(), "Name of vector field")
                    ("fieldz", options::value<std::string>()->required(), "Name of vector field")
                    ("seeds", options::value<long>()->required(), "Number of seed particles")
                    ("steps", options::value<long>()->required(), "Number of Steps")
                    ("length", options::value<float>()->required(), "Length of a single step")
                    ("xrake", options::value<std::string>()->required(), "Number of Steps")
                    ("yrake", options::value<std::string>()->required(), "Number of Steps")
                    ("zrake", options::value<std::string>()->required(), "Number of Steps");
  options::variables_map vm;
  std::ifstream settings_file(std::string(argv[1]), std::ifstream::in);
  options::store(options::parse_config_file(settings_file, desc), vm);
  settings_file.close();
  options::notify(vm);
  if (!(vm.count("datapath")
      && vm.count("fieldx")
      && vm.count("fieldy")
      && vm.count("fieldz")
      && vm.count("seeds")
      && vm.count("steps")
      && vm.count("length")
      && vm.count("xrake")
      && vm.count("yrake")
      && vm.count("zrake")))
  {
    std::cout << "Advection Benchmark" << std::endl << desc << std::endl;
  }

  std::string datapath = vm["datapath"].as<std::string>();
  // Field is not really used
  std::string fieldx = vm["fieldx"].as<std::string>();
  std::string fieldy = vm["fieldy"].as<std::string>();
  std::string fieldz = vm["fieldz"].as<std::string>();

  long numSeeds = vm["seeds"].as<long>();
  long steps = vm["steps"].as<long>();
  float length = vm["length"].as<float>();

  std::cout << "Available # of threads : " << omp_get_max_threads() << std::endl;
  std::cout << "Advection w/ : "
            << "\nData : " << datapath
            << "\nField : " << fieldx << "| " << fieldy << " | " << fieldz
            << "\nSeeds : " << numSeeds
            << "\nSteps : " << steps
            << "\nLength : " << length  << std::endl;

  int res;
  // TODO : read field and variable information
  // Let's assume we have read the data and field
  const std::string filetype = "cf";
  const size_t cache = 2000;
  const size_t threads = 0;

  std::vector<std::string> files;
  GetFilesFromDirectory(datapath, files);
  // no options to set that I know
  std::vector<std::string> fileopts;
  //fileopts.push_back("-project_to_pcs");
  //fileopts.push_back("-vertical_xform");
  VAPoR::DataMgr datamgr(filetype, cache, threads);
  res = datamgr.Initialize(files, fileopts);
  if(res < 0)
  {
    std::cerr << "Failed to intialize CF DataMGR" << std::endl;
    exit(EXIT_FAILURE);
  }
  std::vector<double> timeCoords;
  datamgr.GetTimeCoordinates(timeCoords);
  double initTime = timeCoords.at(0);
  // Create particles
  std::vector<flow::Particle> seeds;
  std::vector<double> xrake;
  multitoken_double(vm["xrake"].as<std::string>(), xrake);
  std::cout << "Rake X : [ " << xrake.at(0) << " --- " << xrake.at(1) << "]" << std::endl;
  std::vector<double> yrake;
  multitoken_double(vm["yrake"].as<std::string>(), yrake);
  std::cout << "Rake Y : [ " << yrake.at(0) << " --- " << yrake.at(1) << "]" << std::endl;
  std::vector<double> zrake;
  multitoken_double(vm["zrake"].as<std::string>(), zrake);
  std::cout << "Rake Z : [ " << zrake.at(0) << " --- " << zrake.at(1) << "]" << std::endl;
  /*vector<double> mind, maxd;
  // Get the extents of the dataset
  // I don't know why this requires a variable name!
  res = datamgr.GetVariableExtents(0, "uinterp", -1, -1, mind, maxd);
  std::cout << "****************" << datamgr.GetNumDimensions("uinterp") << std::endl;
  std::vector<size_t> dims;
  datamgr.GetDimLens("uinterp", dims);
  for(auto& dim : dims)
    std::cout << "**" << dim << std::endl;

  if(res < 0)
  {
    std::cerr << "Failed to retrieve the extents of data" << std::endl;
    exit(EXIT_FAILURE);
  }*/
 // Populate seeds array
  // Random for now, let users configure later
  GenerateSeeds(seeds, xrake, yrake, zrake, numSeeds, initTime);
  std::cout << "Will use " << seeds.size() << " seeds." << std::endl;

  flow::VaporField velocityField(3);
  velocityField.IsSteady = true;
  velocityField.AssignDataManager(&datamgr);
  velocityField.VelocityNames[0] = fieldx.c_str();
  velocityField.VelocityNames[1] = fieldy.c_str();
  velocityField.VelocityNames[2] = fieldz.c_str();

  ParamsBase::StateSave stateSave;// = nullptr;
  VAPoR::FlowParams params(&datamgr, &stateSave);
  params.SetIsSteady(true);
  params.SetSteadyNumOfSteps(steps);
  params.SetFlowDirection(static_cast<int>(VAPoR::FlowDir::FORWARD));
  params.SetSeedGenMode(static_cast<int>(VAPoR::FlowSeedMode::RANDOM));
  params.SetRandomNumOfSeeds(1000);
  velocityField.UpdateParams(&params);

  res = velocityField.CalcDeltaTFromCurrentTimeStep(length);

  // We're building this barebone so no need for field to color with

  // Advect particles in the field
  external::Advection advection;
  advection.UseSeedParticles(seeds);

  auto start = chrono::steady_clock::now();

  int advect = flow::ADVECT_HAPPENED;
  for(size_t step =  advection.GetMaxNumOfPart() - 1;
      step < steps && advect == flow::ADVECT_HAPPENED; steps++)
  {
    advect = advection.AdvectOneStep(&velocityField, length, external::Advection::ADVECTION_METHOD::RK4);
    //std::cout << "Advection happened? " << ((advect != 0) ? "yes" : "no") << std::endl;
  }

  //advect = advection.AdvectSteps(&velocityField, length, steps, external::Advection::ADVECTION_METHOD::RK4);

  auto end = chrono::steady_clock::now();
  const double nanotosec = 1e-9;
  auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - start).count() * nanotosec;
  cout << "Elapsed time : " << elapsed << " sec." << endl;

  PrintStreams(advection);
  /*res = advection.AdvectTillTime(&velocityField, 0, length, 10, flow::Advection::ADVECTION_METHOD::RK4);*/
  //for(auto& seed : seeds)
  //  std::cout << seed.location.x << " : " << seed.time << std::endl;
  return 0;
}
