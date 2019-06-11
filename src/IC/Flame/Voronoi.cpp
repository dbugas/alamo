#include <AMReX_ParallelDescriptor.H>
#include "Voronoi.H"

namespace IC
{
namespace Flame
{
Voronoi::Voronoi (amrex::Vector<amrex::Geometry> &_geom, int _number_of_grains)
	: IC(_geom), number_of_grains(_number_of_grains)
{
amrex::Print() << "number of grains in IC/Flame/Voronoi.cpp  "<<_number_of_grains << "\n";	
AMREX_D_TERM(voronoi_x.resize(number_of_grains);,
					 voronoi_y.resize(number_of_grains);,
					 voronoi_z.resize(number_of_grains););

	for (int n = 0; n<number_of_grains; n++)
		{
			AMREX_D_TERM(voronoi_x[n] = geom[0].ProbLo(0) + (geom[0].ProbHi(0)-geom[0].ProbLo(0))*((amrex::Real)rand())/((amrex::Real)RAND_MAX);,
							 voronoi_y[n] = geom[0].ProbLo(1) + (geom[0].ProbHi(1)-geom[0].ProbLo(1))*((amrex::Real)rand())/((amrex::Real)RAND_MAX);,
							 voronoi_z[n] = geom[0].ProbLo(2) + (geom[0].ProbHi(2)-geom[0].ProbLo(2))*((amrex::Real)rand())/((amrex::Real)RAND_MAX););
		}
}
  
void Voronoi::Add(const int lev, amrex::Vector<amrex::MultiFab * > &field)
{
	AMREX_D_TERM(amrex::Real sizex = geom[0].ProbHi()[0] - geom[0].ProbLo()[0];,
		     amrex::Real sizey = geom[0].ProbHi()[1] - geom[0].ProbLo()[1];,
		     amrex::Real sizez = geom[0].ProbHi()[2] - geom[0].ProbLo()[2];)

		for (amrex::MFIter mfi(*field[lev],true); mfi.isValid(); ++mfi)
			{
				const amrex::Box& box = mfi.tilebox();

				amrex::BaseFab<amrex::Real> &field_box = (*field[lev])[mfi];

				AMREX_D_TERM(for (int i = box.loVect()[0]-field[lev]->nGrow(); i<=box.hiVect()[0]+field[lev]->nGrow(); i++),
								 for (int j = box.loVect()[1]-field[lev]->nGrow(); j<=box.hiVect()[1]+field[lev]->nGrow(); j++),
								 for (int k = box.loVect()[2]-field[lev]->nGrow(); k<=box.hiVect()[2]+field[lev]->nGrow(); k++))
					{
						AMREX_D_TERM(amrex::Real x = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];,
										 amrex::Real y = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];,
										 amrex::Real z = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);

						amrex::IntVect m(AMREX_D_DECL(i,j,k));

						amrex::Real min_distance = std::numeric_limits<amrex::Real>::infinity();
						int min_grain_id = -1;
						field_box(m) = 0.; // initialize
						for (int n = 0; n<number_of_grains; n++)
							{
								amrex::Real d = sqrt(AMREX_D_TERM((x-voronoi_x[n])*(x-voronoi_x[n]), + (y-voronoi_y[n])*(y-voronoi_y[n]), + (z-voronoi_z[n])*(z-voronoi_z[n])));

								if (geom[0].isPeriodic(0))
									{
										d = std::min(d,sqrt(AMREX_D_TERM((x-voronoi_x[n] + sizex)*(x-voronoi_x[n] + sizex), + (y-voronoi_y[n])*(y-voronoi_y[n]), + (z-voronoi_z[n])*(z-voronoi_z[n]))));
										d = std::min(d,sqrt(AMREX_D_TERM((x-voronoi_x[n] - sizex)*(x-voronoi_x[n] - sizex), + (y-voronoi_y[n])*(y-voronoi_y[n]), + (z-voronoi_z[n])*(z-voronoi_z[n]))));
									}
								if (geom[0].isPeriodic(1))
									{
										d = std::min(d,sqrt(AMREX_D_TERM((x-voronoi_x[n])*(x-voronoi_x[n]), + (y-voronoi_y[n] + sizey)*(y-voronoi_y[n] + sizey), + (z-voronoi_z[n])*(z-voronoi_z[n]))));
										d = std::min(d,sqrt(AMREX_D_TERM((x-voronoi_x[n])*(x-voronoi_x[n]), + (y-voronoi_y[n] - sizey)*(y-voronoi_y[n] - sizey), + (z-voronoi_z[n])*(z-voronoi_z[n]))));
									}
#if AMREX_SPACEDIM>2
								if (geom[0].isPeriodic(2))
									{
										d = std::min(d,sqrt(AMREX_D_TERM((x-voronoi_x[n])*(x-voronoi_x[n]), + (y-voronoi_y[n])*(y-voronoi_y[n]), + (z-voronoi_z[n] + sizez)*(z-voronoi_z[n] + sizez))));
										d = std::min(d,sqrt(AMREX_D_TERM((x-voronoi_x[n])*(x-voronoi_x[n]), + (y-voronoi_y[n])*(y-voronoi_y[n]), + (z-voronoi_z[n] - sizez)*(z-voronoi_z[n] - sizez))));
									}
#endif
								if (d<min_distance)
									{
										min_distance = d;
										min_grain_id = n;
									}
							}
						field_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) = (amrex::Real)(min_grain_id);
					}
			}
}
}
}
