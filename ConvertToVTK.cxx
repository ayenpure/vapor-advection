#include <chrono>
#include <filesystem>
#include <iostream>
#include <vector>
#include <random>

#include "boost/program_options.hpp"

#include <vapor/DataMgr.h>
#include <vapor/FileUtils.h>
#include <vapor/VaporField.h>

#include "Advection.h"
#include "FTLEHelper.h"

using namespace VAPoR;


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
  const size_t cache = 2000;
  const size_t threads = 0;

  std::vector<std::string> files;
  {
    namespace fs = std::__fs::filesystem;
    for (const auto & entry : fs::directory_iterator(datapath))
    {
        files.push_back(entry.path());
        std::cout << entry.path() << std::endl;
    }
  }
  printf("file = %s\n", files[0].c_str());
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

  //std::cout << "Number of compression levels : " << datamgr.GetNumRefLevels() << std::endl;

  std::vector<std::string> dimnames = datamgr.GetDimensionNames();
  for(auto& name : dimnames)
  {
    std::cout << " | " << name << "(" << datamgr.IsCompressed(name) << ")" << " | ";
  }
  std::cout << std::endl;

  std::vector<double> timecoords = datamgr.GetTimeCoordinates();
  for(auto& time : timecoords)
    std::cout << " | " << time << " | ";
  std::cout << std::endl;

  // for each file
  {
    // Get coordinate System

    // Get Flow field variables

    // Compose VTK dataset

    // write output with a VTK extension
  }

  return 0;
}
