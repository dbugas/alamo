#include <AMReX_ParallelDescriptor.H>
#include "PackedSpheres.H"

namespace IC
{

namespace Flame
{
PackedSpheres::PackedSpheres (amrex::Vector<amrex::Geometry> &_geom, int _number_of_grains=5)
	: IC(_geom), number_of_grains(5)//(_number_of_grains)
{
//amrex::Print() << "number of grains in IC/Flame/PackedSpheres.cpp  "<<number_of_grains << "\n";
	AMREX_D_TERM(PackedSpheres_x.resize(number_of_grains);,
					 PackedSpheres_y.resize(number_of_grains);,
					 PackedSpheres_z.resize(number_of_grains););
//int a= PackedSpheres_x.resize;
//amrex::Print() << "PackedSpheres_x.resize" <<PackedSpheres_x.resize<< "\n";


	for (int n = 0; n<number_of_grains; n++)
		{
			//AMREX_D_TERM(PackedSpheres_x[n] = geom[0].ProbLo(0) + (geom[0].ProbHi(0)-geom[0].ProbLo(0))*((amrex::Real)rand())/((amrex::Real)RAND_MAX);,
							// PackedSpheres_y[n] = geom[0].ProbLo(1) + (geom[0].ProbHi(1)-geom[0].ProbLo(1))*((amrex::Real)rand())/((amrex::Real)RAND_MAX);,
							// PackedSpheres_z[n] = geom[0].ProbLo(2) + (geom[0].ProbHi(2)-geom[0].ProbLo(2))*((amrex::Real)rand())/((amrex::Real)RAND_MAX););
amrex::Print() << "n " <<n<< "\n";
amrex::Print() << "   " << "\n";

AMREX_D_TERM(PackedSpheres_x[n] = geom[0].ProbLo(0) +((2.5*n)/(15*number_of_grains))*(geom[0].ProbHi(0)-geom[0].ProbLo(0));,
                                                         PackedSpheres_y[n] = geom[0].ProbLo(1) + ((2.5*n)/(5*number_of_grains))*(geom[0].ProbHi(1)-geom[0].ProbLo(1));,
                                                         PackedSpheres_z[n] = geom[0].ProbLo(2) + ((2.5*n)/number_of_grains)*(geom[0].ProbHi(2)-geom[0].ProbLo(2)););


/*amrex::Print() << "geom[0] " <<geom[0]<< "\n";
amrex::Print() << "   " << "\n";*/

/*amrex::Print() << "(amrex::Real)rand()  " <<(amrex::Real)rand()<< "\n";
amrex::Print() << "   " << "\n";*/

/*amrex::Print() << "(amrex::Real)RAND_MAX  " <<(amrex::Real)RAND_MAX<< "\n";
amrex::Print() << "   " << "\n";*/

/*amrex::Print() << "geom[0].ProbLo(0)  " <<geom[0].ProbLo(0)<< "\n";
amrex::Print() << "geom[0].ProbLo(1)  " <<geom[0].ProbLo(1)<< "\n";
amrex::Print() << "geom[0].ProbLo(2)  " <<geom[0].ProbLo(2)<< "\n";
amrex::Print() << "   " << "\n";

amrex::Print() << "geom[0].ProbHi(0)  " <<geom[0].ProbHi(0)<< "\n";
amrex::Print() << "geom[0].ProbHi(1)  " <<geom[0].ProbHi(1)<< "\n";
amrex::Print() << "geom[0].ProbHi(2)  " <<geom[0].ProbHi(2)<< "\n";
amrex::Print() << "   " << "\n";

amrex::Print() << "PackedSpheres_x[n] " <<PackedSpheres_x[n]<< "\n";
amrex::Print() << "PackedSpheres_y[n] " <<PackedSpheres_y[n]<< "\n";
amrex::Print() << "PackedSpheres_z[n] " <<PackedSpheres_z[n]<< "\n";
amrex::Print() << "   " << "\n";*/

	}
}

  
void PackedSpheres::Add(const int lev, amrex::Vector<amrex::MultiFab * > &field)
{
	AMREX_D_TERM(amrex::Real sizex = geom[0].ProbHi()[0] - geom[0].ProbLo()[0];,
		     amrex::Real sizey = geom[0].ProbHi()[1] - geom[0].ProbLo()[1];,
		     amrex::Real sizez = geom[0].ProbHi()[2] - geom[0].ProbLo()[2];)


/*amrex::Print() << "geom[0].ProbLo()[0]  " <<geom[0].ProbLo()[0]<< "\n";
amrex::Print() << "geom[0].ProbLo()[1]  " <<geom[0].ProbLo()[1]<< "\n";
amrex::Print() << "geom[0].ProbLo()[2]  " <<geom[0].ProbLo()[2]<< "\n";

amrex::Print() << "geom[0].ProbHi()[0]  " <<geom[0].ProbHi()[0]<< "\n";
amrex::Print() << "geom[0].ProbHi()[1]  " <<geom[0].ProbHi()[1]<< "\n";
amrex::Print() << "geom[0].ProbHi()[2]  " <<geom[0].ProbHi()[2]<< "\n";

amrex::Print() << "size x " <<sizex<< "\n";
amrex::Print() << "size y " <<sizey<< "\n";
amrex::Print() << "size z " <<sizez<< "\n";*/


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

/*amrex::Print() << "i before    " <<i<< "\n";
amrex::Print() << "   " << "\n";

amrex::Print() << "j before  " <<i<< "\n";
amrex::Print() << "   " << "\n";

amrex::Print() << "k before  " <<i<< "\n";
amrex::Print() << "   " << "\n";

amrex::Print() << "box  " <<box<< "\n";
amrex::Print() << "   " << "\n";

//amrex::Print() << "mfi  " <<mfi<< "\n";
amrex::Print() << "   " << "\n";


amrex::Print() << "lev    " <<lev<< "\n";
amrex::Print() << "   " << "\n";

amrex::Print() << "(*field[lev])[mfi]    " <<(*field[lev])[mfi]<< "\n";
amrex::Print() << "   " << "\n";

amrex::Print() << "box.loVect()[0]    " <<box.loVect()[0]<< "\n";
amrex::Print() << "box.loVect()[1]    " <<box.loVect()[1]<< "\n";
amrex::Print() << "box.loVect()[2]    " <<box.loVect()[2]<< "\n";
amrex::Print() << "   " << "\n";

amrex::Print() << "box.hiVect()[0]    " <<box.hiVect()[0]<< "\n";
amrex::Print() << "box.hiVect()[1]    " <<box.hiVect()[1]<< "\n";
amrex::Print() << "box.hiVect()[2]    " <<box.hiVect()[2]<< "\n";
amrex::Print() << "   " << "\n";


amrex::Print() << "geom[lev].ProbLo()[0]    " <<geom[lev].ProbLo()[0]<< "\n";
amrex::Print() << "geom[lev].ProbLo()[1]    " <<geom[lev].ProbLo()[1]<< "\n";
amrex::Print() << "geom[lev].ProbLo()[2]    " <<geom[lev].ProbLo()[2]<< "\n";
amrex::Print() << "   " << "\n";

amrex::Print() << "(amrex::Real)(i) + 0.5)    " <<(amrex::Real)(i) + 0.5<< "\n";
amrex::Print() << "(amrex::Real)(j) + 0.5   " <<(amrex::Real)(j) + 0.5<< "\n";
amrex::Print() << "(amrex::Real)(k) + 0.5    " <<(amrex::Real)(k) + 0.5<< "\n";
amrex::Print() << "   " << "\n";

amrex::Print() << "i    " <<i<< "\n";
amrex::Print() << "   " << "\n";

amrex::Print() << "j    " <<j<< "\n";
amrex::Print() << "   " << "\n";

amrex::Print() << "k    " <<k<< "\n";
amrex::Print() << "   " << "\n";


amrex::Print() << " geom[lev].CellSize()[0]    " << geom[lev].CellSize()[0]<< "\n";
amrex::Print() << " geom[lev].CellSize()[1]    " << geom[lev].CellSize()[1]<< "\n";
amrex::Print() << " geom[lev].CellSize()[2]    " << geom[lev].CellSize()[2]<< "\n";
amrex::Print() << "   " << "\n"*/

//amrex::Print() << "Real x    " <<x<< "\n";
//amrex::Print() << "Real y    " <<y<< "\n";
//amrex::Print() << "Real z    " <<z<< "\n";
//amrex::Print() << "   " << "\n";
						amrex::IntVect m(AMREX_D_DECL(i,j,k));
//amrex::Print() << "IntVect    " <<m<< "\n";
						amrex::Real min_distance = std::numeric_limits<amrex::Real>::infinity();
//amrex::Print() << "min_distance    " <<min_distance<< "\n";
						int min_grain_id = -1;
						field_box(m) = 0.; // initialize
int theta=0;
int h1=0;
int h2=0;
int r=3;
int i_max=1000;
int inc_theta=(15*3.14)/180;
int max_theta=(360*3.14)/180;




						for (int n = 0; n<number_of_grains; n++)
							{

for(int i=1;i<i_max;i++)
{
while (theta<=max_theta)
{

amrex::Real e1=h1+r*cos(theta);
amrex::Real e2=h2+r*sin(theta);
theta =theta+inc_theta;
i++;
}}
								//amrex::Real d = sqrt(AMREX_D_TERM((x-PackedSpheres_x[n])*(x-PackedSpheres_x[n]), + (y-PackedSpheres_y[n])*(y-PackedSpheres_y[n]), + (z-PackedSpheres_z[n])*(z-PackedSpheres_z[n])));
//amrex::Print() << "PackedSpheres_x[n]    " <<PackedSpheres_x[n] << "\n";
//amrex::Print() << "PackedSpheres_y[n]    " <<PackedSpheres_y[n] << "\n";
//amrex::Print() << "PackedSpheres_z[n]    " <<PackedSpheres_z[n] << "\n";
//amrex::Print() << "Real d    " <<d << "\n";
//amrex::Print() << "field_box(m)    " <<field_box(m) << "\n";

//amrex::Real d = sqrt(AMREX_D_TERM((x-PackedSpheres_x[n] + sizex)*(x-PackedSpheres_x[n] + sizex), + (y-PackedSpheres_y[n]+sizey)*(y-PackedSpheres_y[n]+sizey), + (z-PackedSpheres_z[n]+sizez)*(z-PackedSpheres_z[n]+sizez)));
// d = sqrt(AMREX_D_TERM((x-PackedSpheres_x[n] - sizex)*(x-PackedSpheres_x[n] - sizex), + (y-PackedSpheres_y[n]-sizey)*(y-PackedSpheres_y[n]-sizey), + (z-PackedSpheres_z[n]-sizez)*(z-PackedSpheres_z[n]-sizez)));


/*							if (geom[0].isPeriodic(0))
									{
										d = std::min(d,sqrt(AMREX_D_TERM((x-PackedSpheres_x[n] + sizex)*(x-PackedSpheres_x[n] + sizex), + (y-PackedSpheres_y[n])*(y-PackedSpheres_y[n]), + (z-PackedSpheres_z[n])*(z-PackedSpheres_z[n]))));
//amrex::Print() << "Periodic(0)   d1    " <<d << "\n";		

//amrex::Real d = sqrt(AMREX_D_TERM((x-PackedSpheres_x[n] )*(x-PackedSpheres_x[n] ), + (y-PackedSpheres_y[n])*(y-PackedSpheres_y[n]), + (z-PackedSpheres_z[n])*(z-PackedSpheres_z[n]))));
							
                                                                           	d = std::min(d,sqrt(AMREX_D_TERM((x-PackedSpheres_x[n] - sizex)*(x-PackedSpheres_x[n] - sizex), + (y-PackedSpheres_y[n])*(y-PackedSpheres_y[n]), + (z-PackedSpheres_z[n])*(z-PackedSpheres_z[n]))));

amrex::Print() << "Periodic(0)   d2    " <<d << "\n";

//amrex:: Real d = sqrt(AMREX_D_TERM((x-PackedSpheres_x[n] )*(x-PackedSpheres_x[n]), + (y-PackedSpheres_y[n])*(y-PackedSpheres_y[n]), + (z-PackedSpheres_z[n])*(z-PackedSpheres_z[n])));

					}
								if (geom[0].isPeriodic(1))
									{
										d = std::min(d,sqrt(AMREX_D_TERM((x-PackedSpheres_x[n])*(x-PackedSpheres_x[n]), + (y-PackedSpheres_y[n] + sizey)*(y-PackedSpheres_y[n] + sizey), + (z-PackedSpheres_z[n])*(z-PackedSpheres_z[n]))));
//amrex::Print() << "Periodic(1)   d1    " <<d << "\n";			

//amrex::Real d = sqrt(AMREX_D_TERM((x-PackedSpheres_x[n])*(x-PackedSpheres_x[n]), + (y-PackedSpheres_y[n])*(y-PackedSpheres_y[n]), + (z-PackedSpheres_z[n])*(z-PackedSpheres_z[n])));
	
                                    						d = std::min(d,sqrt(AMREX_D_TERM((x-PackedSpheres_x[n])*(x-PackedSpheres_x[n]), + (y-PackedSpheres_y[n] - sizey)*(y-PackedSpheres_y[n] - sizey), + (z-PackedSpheres_z[n])*(z-PackedSpheres_z[n]))));
								
//amrex::Print() << "Periodic(1)   d2    " <<d << "\n";	

//amrex::Real d = sqrt(AMREX_D_TERM((x-PackedSpheres_x[n])*(x-PackedSpheres_x[n]), + (y-PackedSpheres_y[n])*(y-PackedSpheres_y[n]), + (z-PackedSpheres_z[n])*(z-PackedSpheres_z[n])));
                                      }
#if AMREX_SPACEDIM>2
								if (geom[0].isPeriodic(2))
									{
										d = std::min(d,sqrt(AMREX_D_TERM((x-PackedSpheres_x[n])*(x-PackedSpheres_x[n]), + (y-PackedSpheres_y[n])*(y-PackedSpheres_y[n]), + (z-PackedSpheres_z[n] + sizez)*(z-PackedSpheres_z[n] + sizez))));
//amrex::Print() << "Periodic(2)   d1    " <<d << "\n";			
//amrex::Real d = sqrt(AMREX_D_TERM((x-PackedSpheres_x[n])*(x-PackedSpheres_x[n]), + (y-PackedSpheres_y[n])*(y-PackedSpheres_y[n]), + (z-PackedSpheres_z[n])*(z-PackedSpheres_z[n])));	
                                    						d = std::min(d,sqrt(AMREX_D_TERM((x-PackedSpheres_x[n])*(x-PackedSpheres_x[n]), + (y-PackedSpheres_y[n])*(y-PackedSpheres_y[n]), + (z-PackedSpheres_z[n] - sizez)*(z-PackedSpheres_z[n] - sizez))));
//amrex::Print() << "Periodic(2)   d2    " <<d << "\n";			

//amrex::Real d = sqrt(AMREX_D_TERM((x-PackedSpheres_x[n])*(x-PackedSpheres_x[n]), + (y-PackedSpheres_y[n])*(y-PackedSpheres_y[n]), + (z-PackedSpheres_z[n])*(z-PackedSpheres_z[n])));					
                                	}
#endif*/
//								if (d<min_distance)
//									{
//										min_distance = d;
//										min_grain_id = n;
//									}
							}
//						field_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) = (amrex::Real)(min_grain_id);
					}
			}
}
}
}
