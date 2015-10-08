#include "mfem.hpp"
#include "receivers.hpp"

using namespace mfem;

//==============================================================================
//
// ReceiversSet
//
//==============================================================================
ReceiversSet::ReceiversSet()
  : _n_receivers(0),
    _receivers(),
    _cells_containing_receivers()
{ }

ReceiversSet::ReceiversSet(const ReceiversSet& rec)
  : _n_receivers(rec._n_receivers),
    _receivers(rec._receivers),
    _cells_containing_receivers(rec._cells_containing_receivers)
{ }

ReceiversSet& ReceiversSet::operator =(const ReceiversSet& rec)
{
  _n_receivers = rec._n_receivers;
  _receivers = rec._receivers;
  _cells_containing_receivers = rec._cells_containing_receivers;
  return *this;
}

int find_element(int n_cells_x, int n_cells_y, double sx, double sy,
                 const Vertex &point, bool throw_exception)
{
  const double px = point(0); // coordinates of the point of interest
  const double py = point(1);

  const double x0 = 0.0; // limits of the rectangular mesh
  const double x1 = sx;
  const double y0 = 0.0;
  const double y1 = sy;

  // check that the point is within the mesh
  const double tol = FIND_CELL_TOLERANCE;
  if (px < x0 - tol || px > x1 + tol ||
      py < y0 - tol || py > y1 + tol)
  {
    if (throw_exception)
      MFEM_ABORT("The given point [" + d2s(px) + "," + d2s(py) + "] doesn't "
                 "belong to the rectangular mesh");

    return -1; // to show that the point in not here
  }

  // since the elements of the rectangular mesh are numerated in the following
  // way:
  // -----------
  // | 2  | 3  |
  // -----------
  // | 0  | 1  |
  // -----------
  // we can simplify the search of the element containing the given point:

  // find the x-range where the point is
  int nx = -1;
  const double hx = (x1 - x0) / n_cells_x;
  for (int i = 1; i < n_cells_x+1 && nx == -1; ++i)
    if (px < x0 + i*hx + tol)
      nx = i - 1;
  MFEM_VERIFY(nx != -1, "The x-range is not found");

  // find the y-range where the point is
  int ny = -1;
  const double hy = (y1 - y0) / n_cells_y;
  for (int i = 1; i < n_cells_y+1 && ny == -1; ++i)
    if (py < y0 + i*hy + tol)
      ny = i - 1;
  MFEM_VERIFY(ny != -1, "The y-range is not found");

  return (ny*n_cells_x + nx);
}

void ReceiversSet::
find_cells_containing_receivers(int nx, int ny, double sx, double sy)
{
  MFEM_VERIFY(!_receivers.empty(), "The receivers hasn't been distributed yet");

  _cells_containing_receivers.clear();
  _cells_containing_receivers.resize(_n_receivers);

  // we throw an exception if we don't find a cell containing a receiver
  const bool throw_exception = true;

  for (int p = 0; p < _n_receivers; ++p)
    _cells_containing_receivers[p] = find_element(nx, ny, sx, sy, _receivers[p],
                                                  throw_exception);
}



//==============================================================================
//
// ReceiversLine
//
//==============================================================================
ReceiversLine::ReceiversLine()
  : _start(),
    _end()
{ }

ReceiversLine::ReceiversLine(const ReceiversLine& rec)
  : ReceiversSet(rec),
    _start(rec._start),
    _end(rec._end)
{ }

ReceiversLine& ReceiversLine::operator =(const ReceiversLine& rec)
{
  ReceiversSet::operator =(rec);
  _start = rec._start;
  _end   = rec._end;
  return *this;
}

void ReceiversLine::init(std::ifstream &in)
{
  MFEM_VERIFY(in.is_open(), "The stream for reading receivers is not open");

  std::string tmp;

  in >> _n_receivers; getline(in, tmp);
  in >> _start(0); getline(in, tmp);
  in >> _end(0);   getline(in, tmp);
  in >> _start(1); getline(in, tmp);
  in >> _end(1);   getline(in, tmp);

  MFEM_VERIFY(_n_receivers > 0, "The number of receivers (" + d2s(_n_receivers)+
              ") must be >0");
}

void ReceiversLine::distribute_receivers()
{
  const double x0 = _start(0);
  const double x1 = _end(0);
  const double y0 = _start(1);
  const double y1 = _end(1);

  _receivers.resize(_n_receivers);

  const double dx = (x1 - x0) / (_n_receivers-1);
  const double dy = (y1 - y0) / (_n_receivers-1);

  for (int i = 0; i < _n_receivers; ++i)
  {
    const double x = (i == _n_receivers-1 ? x1 : x0 + i*dx);
    const double y = (i == _n_receivers-1 ? y1 : y0 + i*dy);
    _receivers[i] = Vertex(x, y);
  }
}

std::string ReceiversLine::description() const
{
  return "_rec_line_x" + d2s(_start(0)) + "_" + d2s(_end(0)) +
         "_y" + d2s(_start(1)) + "_" + d2s(_end(1));
}



