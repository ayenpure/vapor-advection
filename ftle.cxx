#include <chrono>
#include <iostream>
#include <vector>
#include <random>

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

#include <vapor/DataMgr.h>
#include <vapor/FileUtils.h>
#include <vapor/VaporField.h>

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

/*void GenerateSeeds(std::vector<flow::Particle>& seeds,
                   const std::vector<double>& rake,
                   const long numOfSeeds)
{
  // We need uniform seeding here throughout the grid
  // Assume we Have that
  std::cout << "[begin] generate seeds" << std::endl;
  std::cout << "Rake size : " << rake.size() << std::endl;
  const unsigned int randSeed = 32;
  std::mt19937 gen(randSeed); //Standard mersenne_twister_engine
  std::uniform_real_distribution<double> distX( rake.at(0), rake.at(1));
  std::uniform_real_distribution<double> distY( rake.at(2), rake.at(3));
  std::uniform_real_distribution<double> distZ( rake.at(4), rake.at(5));
  const float timeVal = 0.;
  seeds.resize(numOfSeeds);
  for( long i = 0; i < numOfSeeds; i++ )
  {
    seeds[i].location.x = distX(gen);
    seeds[i].location.y = distY(gen);
    seeds[i].location.z = distZ(gen);
    seeds[i].time       = timeVal;
  }
  std::cout << "[end] generate seeds" << std::endl;
}*/

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
                    ("seeds", options::value<long>()->required(), "Number of seed particles")
                    ("steps", options::value<long>()->required(), "Number of Steps")
                    ("length", options::value<float>()->required(), "Length of a single step");
  options::variables_map vm;
  options::store(options::parse_command_line(argc, argv, desc), vm); // can throw
  options::notify(vm);
  if (!(vm.count("data")
      && vm.count("field")
      && vm.count("seeds")
      && vm.count("steps")
      && vm.count("length")))
  {
    std::cout << "Advection Benchmark" << std::endl << desc << std::endl;
  }

  std::string datapath = vm["data"].as<std::string>();
  // Field is not really used
  std::string field = vm["field"].as<std::string>();
  long numSeeds = vm["seeds"].as<long>();
  long steps = vm["steps"].as<long>();
  float length = vm["length"].as<float>();

  std::cout << "Advection w/ : "
            << "\nData : " << datapath
            << "\nField : " << field
            << "\nSteps : " << steps
            << "\nLength : " << length  << std::endl;

  int res;
  // TODO : read field and variable information
  // Let's assume we have read the data and field
  const std::string filetype = "cf";
  const size_t cache = 10000;
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

  std::vector<double> timecoords;
  datamgr.GetTimeCoordinates(timecoords);
  for(auto& time : timecoords)
    std::cout << time << ", ";
  std::cout << std::endl;

  // Create particles
  std::vector<flow::Particle> seeds;
  std::vector<double> rake;
  vector<double> mind, maxd;
  // Get the extents of the dataset
  // I don't know why this requires a variable name!
  res = datamgr.GetVariableExtents(0, "u", -1, -1, mind, maxd);
  if(res < 0)
  {
    std::cerr << "Failed to retrieve the extents of data" << std::endl;
    exit(EXIT_FAILURE);
  }
  //Build a seeding rake 
  rake.push_back(mind.at(0));
  rake.push_back(maxd.at(0));
  rake.push_back(mind.at(1));
  rake.push_back(maxd.at(1));
  rake.push_back(mind.at(2));
  rake.push_back(maxd.at(2));

  std::vector<size_t> dims;
  datamgr.GetDimLens("u", dims);

  detail::GridMetaData metaData(dims, rake);
  std::cout << "dims" << std::endl;
  for(auto& dim : dims)
    std::cout << dim << ", ";
  std::cout << std::endl;
  std::cout << "rake" << std::endl;
  for(auto& _rake : rake)
    std::cout << _rake << ", ";
  std::cout << std::endl;
  // Populate seeds array
  metaData.GetSeeds(seeds);
  std::cout << "Will use " << seeds.size() << " seeds." << std::endl;

  flow::VaporField velocityField(8);
  velocityField.IsSteady = true;
  velocityField.AssignDataManager(&datamgr);
  velocityField.VelocityNames[0] = "u";
  velocityField.VelocityNames[1] = "v";
  velocityField.VelocityNames[2] = "w";

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
  //res = advection.AdvectTillTime(&velocityField, 0, length, 10, external::Advection::ADVECTION_METHOD::RK4);
  // Extract streams from the advection class
  std::vector<flow::Particle> endLocations;
  size_t streams = advection.GetNumberOfStreams();
  for(size_t index = 0; index < streams; index++)
  {
    endLocations.push_back(advection.GetStreamAt(index).back());
  }
  /*We'll get streams as an output over here*/
  /*Future optimization : if only FTLE is required, do not calculate the streams*/
  /*FTLE steps*/
  // 1. Calculate Gradiant
  // 2. Calculate Caunchy Green Tensor
  // 3. Calculate exponents
  // 4. Allow rich set of exponents calculation
  std::vector<double> FTLEfield;
  CalculateFTLE(seeds, endLocations, metaData, 10, FTLEfield);

  auto end = chrono::steady_clock::now();
  const double nanotosec = 1e-9;
  auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - start).count() * nanotosec;
  cout << "Elapsed time in nanoseconds : " << elapsed << " sec." << endl;

  return 0;
}
