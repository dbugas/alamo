#include <iostream>
#include <fstream>
#include <iomanip>

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_Array.H>

#include "Util/Util.H"
#include "Integrator/Flame.H"

int main (int argc, char* argv[])
{

	Util::Initialize(argc,argv);
	amrex::Print() << "Hello World n flame.cc"<< "\n";
	Util::Message(INFO,"Hello World in flame.cc");
	Integrator::Integrator *flame =
		new Integrator::Flame();
	flame->InitData();
	flame->Evolve();
	delete flame;

	Util::Finalize();
} 
