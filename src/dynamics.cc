#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "Integrator/ElasticDynamics.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"

int main (int argc, char* argv[])
{
	Util::Initialize(argc,argv);

	srand(1);
	Integrator::Integrator *model =
		new Integrator::ElasticDynamics();
	model->InitData();
	model->Evolve();
	delete model;

	Util::Finalize();
} 
