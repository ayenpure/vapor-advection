#ifndef GRID_META_DATA
#define GRID_META_DATA

#include <vector>

namespace detail
{

class GridMetaData
{
public:
  GridMetaData(const std::vector<long>& _dims)
  : dims(_dims)
  {
    this->PlaneSize = dims.at(0)*dims.at(1);
    this->RowSize = dims.at(0);

    if(this->dims[2] == 1)
      _2DGrid = true;
  }

  long long int GetNumberOfPoints() const
  {
    std::cout <<  dims[0] * dims[1] * dims[2] << std::endl;
    return dims[0] * dims[1] * dims[2];
  }

  void GetSeeds(std::vector<flow::Particle>& seeds) const
  {
    // Cartesian product of x, y, z coords
    size_t xind, yind, zind;
    long long int index = 0;
    const double timeVal = 0.;
    seeds.resize(this->GetNumberOfPoints());
    // Z grows slowest
    for(zind = 0; zind < dims.at(2); zind++)
    {
      for(yind = 0; yind < dims.at(1); yind++)
      {
        // X grows fastest
        for(xind = 0; xind < dims.at(0); xind++)
        {
          seeds[index].location.x = xCoords.at(xind);
          seeds[index].location.y = yCoords.at(yind);
          seeds[index].location.z = zCoords.at(zind);
          seeds[index].time       = timeVal;
          ++index;
        }
      }
    }
  }

  void GetLogicalIndex(const long long int index, long* logical) const
  {
    logical[0] = index % dims[0];
    logical[1] = (index / dims[0]) % dims[1];
    if(this->_2DGrid)
      return;
    logical[2] = index / (dims[0] * dims[1]);
  }

  void GetNeighborIndices(const long long int index, long long int*  neighbors) const
  {
    long logical[3];
    GetLogicalIndex(index, logical);

    // For differentials w.r.t delta in x
    neighbors[0] = (logical[0] == 0) ? index : index - 1;
    neighbors[1] = (logical[0] == dims[0] - 1) ? index : index + 1;
    // For differentials w.r.t delta in y
    neighbors[2] = (logical[1] == 0) ? index : index - RowSize;
    neighbors[3] = (logical[1] == dims[1] - 1) ? index : index + RowSize;
    if(this->_2DGrid)
      return;
    // For differentials w.r.t delta in z
    neighbors[4] = (logical[2] == 0) ? index : index - PlaneSize;
    neighbors[5] = (logical[2] == dims[2] - 1) ? index : index + PlaneSize;
  }

  bool IsTwoDimentional() const
  {
    return this->_2DGrid;
  }

private:
  std::vector<long> dims;
  long long int PlaneSize;
  long long int RowSize;
  bool _2DGrid = false;
};

}

#endif
