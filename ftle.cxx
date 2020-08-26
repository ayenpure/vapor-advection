#include <chrono>
#include <iostream>
#include <limits>
#include <vector>
#include <random>

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

#include <vapor/Advection.h>
#include <vapor/DataMgr.h>
#include <vapor/FileUtils.h>
#include <vapor/VaporField.h>

#include "GridMetaData.h"
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

void GenerateSeeds(std::vector<flow::Particle>& seeds,
                   const std::vector<double>& overlapmin,
                   const std::vector<double>& overlapmax,
                   const std::vector<long>& dimensions,
                   const double initTime)
{
    /* Create arrays that contain X, Y, and Z coordinates */
    float start[3], step[3];
    for( int i = 0; i < dimensions.size(); i++ )    // for each of the X, Y, Z dimensions
    {
        if( dimensions[i] == 1 )     // one seed in this dimension
        {
            start[i] = overlapmin[i] + 
                       0.5f * (overlapmax[i] - overlapmin[i]);
            step[i]  = 0.0f;
        }
        else                        // more than one seed in this dimension
        {
            start[i] = overlapmin[i];
            step[i]  = (overlapmax[i] - overlapmin[i]) / float(dimensions[i] - 1);
        }
    }
    if( dimensions.size() == 2 || dimensions.at(2) == 1)  // put default Z values
    {
        start[2] = 0.0f;
        step[2]  = 0.0f;
    }

    /* Populate the list of seeds */
    float timeVal = initTime;  // Default time value
    glm::vec3 loc;
    seeds.clear();
    long seedsZ;
    if( dimensions.size() == 2 || dimensions.at(2) == 1) seedsZ = 1;
    else            seedsZ = dimensions[2];
    // Reserve enough space at the beginning for performance considerations
    seeds.reserve( seedsZ * dimensions[1] * dimensions[0] );
    for( long k = 0; k < seedsZ; k++ )
        for( long j = 0; j < dimensions[1]; j++ )
            for( long i = 0; i < dimensions[0]; i++ )
            {
                loc.x = start[0] + float(i) * step[0];
                loc.y = start[1] + float(j) * step[1];
                loc.z = start[2] + float(k) * step[2];
                seeds.emplace_back( loc, timeVal );
            }

}

int PrintStreams(flow::Advection& advector)
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

void AdjustOverlap(std::vector<double>& overlapmin,
                   std::vector<double>& overlapmax,
                   const std::vector<double>& mind,
                   const std::vector<double>& maxd)
{
  overlapmin[0] = std::max(overlapmin[0], mind[0]);
  overlapmin[1] = std::max(overlapmin[1], mind[1]);
  overlapmin[2] = std::max(overlapmin[2], mind[2]);

  overlapmax[0] = std::min(overlapmax[0], maxd[0]);
  overlapmax[1] = std::min(overlapmax[1], maxd[1]);
  overlapmax[2] = std::min(overlapmax[2], maxd[2]);
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
                    ("steps", options::value<long>()->required(), "Number of steps")
                    ("length", options::value<float>()->required(), "Length of a single step")
                    ("duration", options::value<double>()->required(), "Duration for advection")
                    ("dimx", options::value<long>()->required(), "Number of seeds in X dimension")
                    ("dimy", options::value<long>()->required(), "Number of seeds in Y dimension")
                    ("dimz", options::value<long>()->required(), "Number of seeds in Z dimension");
  options::variables_map vm;
  std::ifstream settings_file(std::string(argv[1]), std::ifstream::in);
  options::store(options::parse_config_file(settings_file, desc), vm);
  settings_file.close();
  options::notify(vm);
  if (!(vm.count("datapath")
      && vm.count("fieldx")
      && vm.count("fieldy")
      && vm.count("fieldz")
      && vm.count("steps")
      && vm.count("length")
      && vm.count("duration")
      && vm.count("dimx")
      && vm.count("dimy")
      && vm.count("dimz")))
  {
    std::cout << "Advection Benchmark" << std::endl << desc << std::endl;
  }

  const std::string datapath = vm["datapath"].as<std::string>();
  // Field is not really used
  const std::string fieldx = vm["fieldx"].as<std::string>();
  const std::string fieldy = vm["fieldy"].as<std::string>();
  const std::string fieldz = vm["fieldz"].as<std::string>();
  const long steps = vm["steps"].as<long>();
  float length = vm["length"].as<float>();
  const double duration = vm["duration"].as<double>();
  std::vector<long> dimensions;
  const long dimx = vm["dimx"].as<long>();
  dimensions.push_back(dimx);
  const long dimy = vm["dimy"].as<long>();
  dimensions.push_back(dimy);
  const long dimz = vm["dimz"].as<long>();
  dimensions.push_back(dimz);

  std::cout << "Advection w/ : "
            << "\nData : " << datapath
            << "\nField : " << fieldx << "| " << fieldy << " | " << fieldz
            << "\nSteps : " << steps
            << "\nLength : " << length
            << "\nField : " << dimx << "| " << dimy << " | " << dimz << std::endl;

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
  std::cout << "Total time length : "
            << (timeCoords.at(timeCoords.size() - 1) - timeCoords.at(0)) << std::endl;

  // Create particles
  std::vector<flow::Particle> seeds;
  // Get the extents of the dataset
  // I don't know why this requires a variable name!
  vector<double> mind, maxd;
  vector<double> overlapmin(3, std::numeric_limits<double>::min());
  vector<double> overlapmax(3, std::numeric_limits<double>::max());
  res = datamgr.GetVariableExtents(0, fieldx, 0, 0, mind, maxd);
  AdjustOverlap(overlapmin, overlapmax, mind, maxd);
  res = datamgr.GetVariableExtents(0, fieldy, 0, 0, mind, maxd);
  AdjustOverlap(overlapmin, overlapmax, mind, maxd);
  res = datamgr.GetVariableExtents(0, fieldz, 0, 0, mind, maxd);
  AdjustOverlap(overlapmin, overlapmax, mind, maxd);
  std::cout << "Brick Origin : " << overlapmin.at(0) << " " <<overlapmin.at(1) << " " << overlapmin.at(2) << std::endl;
  std::cout << "Brick Size : " << overlapmax.at(0) << " " <<overlapmax.at(1) << " " << overlapmax.at(2) << std::endl;
  // Populate seeds array
  // Random for now, let users configure later
  GenerateSeeds(seeds, overlapmin, overlapmax, dimensions, initTime);
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
  params.SetRefinementLevel(0);
  params.SetCompressionLevel(0);
  params.GetBox()->SetExtents(overlapmin, overlapmax);
  velocityField.UpdateParams(&params);

  res = velocityField.CalcDeltaTFromCurrentTimeStep(length);

  // Advect particles in the field
  flow::Advection advection;
  advection.UseSeedParticles(seeds);

  int advect = flow::ADVECT_HAPPENED;
  advect = advection.AdvectTillTime(&velocityField, initTime, length, (initTime + duration), flow::Advection::ADVECTION_METHOD::RK4);

  auto start = chrono::steady_clock::now();

  // Extract streams from the advection class
  std::vector<flow::Particle> endLocations;
  size_t streams = advection.GetNumberOfStreams();
  for(size_t index = 0; index < streams; index++)
  {
    endLocations.push_back(advection.GetStreamAt(index).back());
  }

  /*We'll get streams as an output over here*/
  /*Future optimization : if only FTLE is required, do not calculate the streams*/
  std::vector<double> bounds;
  bounds.push_back(overlapmin.at(0));
  bounds.push_back(overlapmax.at(0));
  bounds.push_back(overlapmin.at(1));
  bounds.push_back(overlapmax.at(1));
  bounds.push_back(overlapmin.at(2));
  bounds.push_back(overlapmax.at(2));
  detail::GridMetaData metaData(dimensions, bounds);
  std::vector<double> FTLEfield;
  CalculateFTLE(seeds, endLocations, metaData, duration, FTLEfield);

  auto end = chrono::steady_clock::now();
  const double nanotosec = 1e-9;
  auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - start).count() * nanotosec;
  cout << "Elapsed time : " << elapsed << " sec." << endl;

  ofstream fout;
  fout.open("FTLEoutput.dat", ios::binary);
  fout.write(reinterpret_cast<const char*>(&FTLEfield[0]), FTLEfield.size()*sizeof(double));
  fout.close();

  PrintStreams(advection);
  // *res = advection.AdvectTillTime(&velocityField, 0, length, 10, flow::Advection::ADVECTION_METHOD::RK4);*/
  // for(auto& seed : seeds)
  //   std::cout << seed.location.x << " : " << seed.time << std::endl;
  return 0;
}
