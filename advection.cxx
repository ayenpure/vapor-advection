#include <vector>

#include <vapor/Advection.h>
#include <vapor/DataMgr.h>
#include <vapor/FileUtils.h>
#include <vapor/OptionParser.h>
#include <vapor/VaporField.h>

using namespace Wasp;
using namespace VAPoR;

struct
{
  std::string datafile;
  std::string fieldname;
  long numberSteps;
  double stepLength;
} opt;

OptionParser::OptDescRec_T set_opts[] = {
 {"file" , 1, "filename.nc", "Data file : NetCDF?"},
 {"field" , 1, "(x,y,z)?", "Name of the velocity field"},
 {"steps", 1, "100", "Number of maximum number of steps"},
 {"length", 1, "0.01", "Length of a single step (deltaT)"},
 {nullptr}
};

OptionParser::Option_T get_options[] = {
  {"file", Wasp::CvtToCPPStr, &opt.datafile, sizeof(opt.datafile)},
  {"field", Wasp::CvtToCPPStr, &opt.fieldname, sizeof(opt.fieldname)},
  {"steps", Wasp::CvtToInt, &opt.numberSteps, sizeof(opt.numberSteps)},
  {"length", Wasp::CvtToDouble, &opt.stepLength, sizeof(opt.stepLength)},
  {nullptr}
};

void InitializeOptions (int argc, char **argv, OptionParser options)
{
  std::string ProgName = FileUtils::LegacyBasename(argv[0]);
  if (options.AppendOptions(set_opts) < 0)
  {
    std::cerr << ProgName << " : " << options.GetErrMsg();
    exit(EXIT_FAILURE);
  }
  if (options.ParseOptions(&argc, argv, get_options) < 0)
  {
    std::cerr << ProgName << " : " << options.GetErrMsg();
    exit(EXIT_FAILURE);
  }
}

int main (int argc, char** argv)
{
  OptionParser argopts;
  InitializeOptions(argc, argv, argopts);
  // Input params
  // 1. File
  // 2. Field variable
  // 3. Number of steps
  // 4. Step length (deltaT) : step length is calculated dynamically
  //    This is suggestive?
  // 5. Seeding parameters

  // TODO : read field and variable information
  // Let's assume we have read the data and field
  const std::string fileType = "cf";
  const size_t cache = 200;
  const size_t threads = 2;

  std::vector<std::string> files;
  files.push_back(opt.datafile);

  // no options to set that I knowi
  std::vector<std::string> options;
  VAPoR::DataMgr datamgr("cf", cache, threads);
  int res = datamgr.Initialize(files, options);
  if(res < 0)
  {
    std::cerr << "Failed to intialize CF DataMGR" << std::endl;
    exit(EXIT_FAILURE);
  }
  //VAPoR::FlowParams params;

  std::vector<std::string> variables = datamgr.GetDataVarNames();
  
  for(std::vector<std::string>::iterator it = variables.begin() ; it != variables.end(); ++it)
    std::cout << *it << ',';
  std::cout << '\n';

/*  //Build a seeding rake

  // Create particles
  std::vector<flow::Particle> seeds;
  // Populate seeds array;

  flow::Advection advection;

  flow::VaporField velocityField;
  velocityField.AssignDataManager(&datamgr);
  velocityField.UpdateParams(&params);

  int res;
  res = velocityField.CalcDeltaTFromCurrentTimeStep(deltaT);
*/
  // We're building this barebone so no need for field to color with

  // Advect particles in the field

  return 0;
}
