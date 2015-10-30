#include "elastic_wave2D.hpp"
#include "GLL_quadrature.hpp"
#include "parameters.hpp"
#include "receivers.hpp"
#include "utilities.hpp"

#include <fstream>

using namespace std;
using namespace mfem;



void ElasticWave2D::run()
{
  if (!strcmp(param.method, "fem") || !strcmp(param.method, "FEM"))
    run_FEM_ALID();
  else if (!strcmp(param.method, "sem") || !strcmp(param.method, "SEM"))
    run_SEM_SRM();
  else MFEM_ABORT("Unknown method to be used: " + string(param.method));
}



//------------------------------------------------------------------------------
//
// Auxiliary useful functions
//
//------------------------------------------------------------------------------
Vector compute_function_at_point(double sx, double sy, int nx, int ny,
                                 const Mesh& mesh, const Vertex& point,
                                 int cell, const GridFunction& U)
{
  const double hx = sx / nx;
  const double hy = sy / ny;

  const Element *element = mesh.GetElement(cell);
  Array<int> vert_indices;
  element->GetVertices(vert_indices);
  const double *vert0 = mesh.GetVertex(vert_indices[0]); // min coords

  IntegrationPoint ip;
  ip.x = (point(0) - vert0[0]) / hx;
  ip.y = (point(1) - vert0[1]) / hy;

  Vector values;
  U.GetVectorValue(cell, ip, values);

  return values;
}



Vector compute_function_at_points(double sx, double sy, int nx, int ny,
                                  const Mesh& mesh,
                                  const vector<Vertex>& points,
                                  const vector<int>& cells_containing_points,
                                  const GridFunction& U)
{
  MFEM_ASSERT(points.size() == cells_containing_points.size(), "Sizes mismatch");
  Vector U_at_points(2*points.size());

  for (size_t p = 0; p < points.size(); ++p)
  {
    Vector values = compute_function_at_point(sx, sy, nx, ny, mesh, points[p],
                                              cells_containing_points[p], U);
    MFEM_ASSERT(values.Size() == 2, "Unexpected vector size");
    U_at_points(2*p+0) = values(0);
    U_at_points(2*p+1) = values(1);
  }
  return U_at_points;
}



void output_snapshots(int time_step, const string& snapshot_filebase,
                      const Parameters& param, const GridFunction& U,
                      const GridFunction& V)
{
  Vector u_x, u_y, v_x, v_y;
  U.GetNodalValues(u_x, 1); // components of displacement
  U.GetNodalValues(u_y, 2);
  V.GetNodalValues(v_x, 1); // components of velocity
  V.GetNodalValues(v_y, 2);

  string tstep = d2s(time_step,0,0,0,6), fname;
  if (!strcmp(param.snapshot_format, "bin"))
  {
    fname = snapshot_filebase + "_Ux_t" + tstep + ".bin";
    write_binary(fname.c_str(), u_x.Size(), u_x);
    fname = snapshot_filebase + "_Uy_t" + tstep + ".bin";
    write_binary(fname.c_str(), u_y.Size(), u_y);
    fname = snapshot_filebase + "_Vx_t" + tstep + ".bin";
    write_binary(fname.c_str(), v_x.Size(), v_x);
    fname = snapshot_filebase + "_Vy_t" + tstep + ".bin";
    write_binary(fname.c_str(), v_y.Size(), v_y);
  }
  else if (!strcmp(param.snapshot_format, "vts"))
  {
    fname = snapshot_filebase + "_U_t" + tstep + ".vts";
    write_vts_vector(fname, "U", param.grid.sx, param.grid.sy, param.grid.nx,
                     param.grid.ny, u_x, u_y);
    fname = snapshot_filebase + "_V_t" + tstep + ".vts";
    write_vts_vector(fname, "V", param.grid.sx, param.grid.sy, param.grid.nx,
                     param.grid.ny, v_x, v_y);
  }
  else MFEM_ABORT("Unknown snapshot format: " + string(param.snapshot_format));
}



void output_seismograms(const Parameters& param, const Mesh& mesh,
                        const GridFunction &U, const GridFunction &V,
                        ofstream *seisU, ofstream *seisV)
{
  // for each set of receivers
  for (size_t rec = 0; rec < param.sets_of_receivers.size(); ++rec)
  {
    for (int c = 0; c < N_ELAST_COMPONENTS; ++c)
    {
      MFEM_ASSERT(seisU[rec*N_ELAST_COMPONENTS+c].is_open(), "The stream for "
                  "writing displacement seismograms is not open");
      MFEM_ASSERT(seisV[rec*N_ELAST_COMPONENTS+c].is_open(), "The stream for "
                  "writing velocity seismograms is not open");
    }

    const ReceiversSet *rec_set = param.sets_of_receivers[rec];
    // displacement at the receivers
    const Vector u = compute_function_at_points(param.grid.sx, param.grid.sy,
                                                param.grid.nx, param.grid.ny, mesh,
                                                rec_set->get_receivers(),
                                                rec_set->get_cells_containing_receivers(),
                                                U);
    // velocity at the receivers
    const Vector v = compute_function_at_points(param.grid.sx, param.grid.sy,
                                                param.grid.nx, param.grid.ny, mesh,
                                                rec_set->get_receivers(),
                                                rec_set->get_cells_containing_receivers(),
                                                V);

    MFEM_ASSERT(u.Size() == N_ELAST_COMPONENTS*rec_set->n_receivers() &&
                u.Size() == v.Size(), "Sizes mismatch");

    float val;
    for (int i = 0; i < u.Size(); i += N_ELAST_COMPONENTS)
    {
      for (int j = 0; j < N_ELAST_COMPONENTS; ++j)
      {
        val = u(i+j); // displacement
        seisU[rec*N_ELAST_COMPONENTS + j].write(reinterpret_cast<char*>(&val),
                                                sizeof(val));

        val = v(i+j); // velocity
        seisV[rec*N_ELAST_COMPONENTS + j].write(reinterpret_cast<char*>(&val),
                                                sizeof(val));
      }
    }
  }
}


