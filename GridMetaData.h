#include <vector>

namespace detail
{

class GridMetaData
{
public:
  GridMetaData(long* _dims, std::vector<double>& _bounds)
  , bounds(_bounds)
  {
    this->dims[0] = _dims[0];
    this->dims[1] = _dims[1];
    this->dims[2] = _dims[2];

    this->PlaneSize = _dims[0]*_dims[1];
    this->RowSize = _dims[0];

    // Populate x, y, z coords
    double xSpace = (bounds.at(1) - bounds.at(0)) / (dims[0] - 1);
    double ySpace = (bounds.at(3) - bounds.at(2)) / (dims[1] - 1);
    double zSpace = (bounds.at(5) - bounds.at(4)) / (dims[2] - 1);
    double init = bounds.at(0);
    xCoords.resize(dims[0]);
    yCoords.resize(dims[1]);
    zCoords.resize(dims[2]);
    size_t index;
    for(index = 0; index < dims[0]; index++)
      xCoords[index] = bounds.at(0) + index * xSpace;
    for(index = 0; index < dims[1]; index++)
      yCoords[index] = bounds.at(2) + index * ySpace;
    for(index = 0; index < dims[2]; index++)
      zCoords[index] = bounds.at(4) + index * zSpace;
  }

  void GetSeeds(std::vector<flow::Particle>& seeds) const
  {
    // Cartesian product of x, y, z coords
    size_t xind, yind, zind, index;
    const double timeVal = 0.;
    seeds.resize(this->GetNumberOfPoints());
    // Z grows slowest
    for(zind = 0; zind < dims[2]; zind++)
    {
      for(yind = 0; yind < dims[1]; yind++)
      {
        // X grows fastest
        for(xind = 0; xind < dims[0]; xind++)
        {
          seeds[ind].location.x = xCoords.at(xind);
          seeds[ind].location.y = yCoords.at(yind);
          seeds[ind].location.z = zCoords.at(zind);
          seeds[ind].time       = timeVal;
          ++ind;
        }
      }
    }
  }

  void GetLogicalIndex(const long long int index, long* logical) const
  {
    logical[0] = index % dims[0];
    logical[1] = (index / dims[0]) % dims[1];
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
    // For differentials w.r.t delta in z
    neighbors[4] = (logical[2] == 0) ? index : index - PlaneSize;
    neighbors[5] = (logical[2] == dims[2] - 1) ? index : index + PlaneSize;
  }

  void GetNumberOfPoints() const
  {
    return _dims[0] * _dims[1] * _dims[2];
  }

private:
  long dims[3];
  long long int PlaneSize;
  long long int RowSize;
  std::vector<double> bounds;
  std::vector<double> xCoords;
  std::vector<double> yCoords;
  std::vector<double> zCoords;
};

}
