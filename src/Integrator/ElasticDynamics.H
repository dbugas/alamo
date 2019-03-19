///
/// \file PhaseFieldMicrostructure.H
///
#ifndef INTEGRATOR_ELASTICDYNAMICS_H
#define INTEGRATOR_ELASTICDYNAMICS_H

#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"
#include <AMReX_MLMG.H>

#include "Integrator/Integrator.H"

#include "BC/BC.H"
#include "BC/Constant.H"
#include "IC/TabulatedInterface.H"
#include "IC/PerturbedInterface.H"
#include "IC/Voronoi/Voronoi.H"
#include "IC/Circle.H"

#include "Model/Interface/GrainBoundary/GrainBoundary.H"
#include "Model/Interface/GrainBoundary/Sin.H"
#include "Model/Interface/GrainBoundary/AbsSin.H"
#include "Model/Interface/GrainBoundary/Read.H"

#include "Operator/CellElastic/PolyCrystal/PolyCrystal.H"
#include "Operator/CellElastic/Isotropic/Isotropic.H"
#include "Operator/CellElastic/PolyCrystal/Isotropic/Isotropic.H"
#include "Operator/CellElastic/PolyCrystal/Cubic/Cubic.H"
#include "Operator/CellElastic/Cubic/Cubic.H"
#include "Model/Solid/LinearElastic/Cubic.H"
#include "Model/Solid/LinearElastic/Isotropic.H"
#include "Model/Solid/LinearElastic/Laplacian.H"
#include "Operator/Elastic.H"

namespace Integrator
{

using model_type = Model::Solid::LinearElastic::Isotropic;


///
/// \class PhaseFieldMicrostructure
/// \brief Microstructure evolution with grain boundary anisotropy
///
/// Solve the Allen-Cahn evolution equation for microstructure with parameters \f$\eta_1\ldots\eta_n\f$,
/// where n corresponds to the number of grains.
///
class ElasticDynamics : public Integrator
{
public:
	ElasticDynamics();

protected:

	/// \fn    Advance
	/// \brief Evolve phase field in time
	void Advance (int lev, Real time, Real dt) override;

	void Initialize (int lev) override;

	void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

	void TimeStepBegin(amrex::Real time, int iter) override;
	void TimeStepComplete(amrex::Real time, int iter) override;
	//void Integrate(int amrlev, Set::Scalar time, int step,
	//	       const amrex::MFIter &mfi, const amrex::Box &box) override;

private:

	int number_of_ghost_cells = 2;

	amrex::Vector<std::unique_ptr<amrex::MultiFab> > displacement; ///< Multicomponent field variable storing pointwise displacements
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > displacement_nm1; ///< Multicomponent field variable storing pointwise displacement for t=t_{n-1}
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > displacement_nm2; ///< Multicomponent field variable storing pointwise displacement for t=t_{n-2}
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > velocity; ///< Multicomponent field variable storing pointwise velocity for t=t_n
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > residual; ///< Multicomponent field variable storing pointwise displacements
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > strain; ///< Multicomponent field variable storing pointwise strains (6 components)
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > stress; ///< Multicomponent field variable storing pointwise strains (6 components)
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > divsigma; ///< Multicomponent field variable storing pointwise strains (6 components)
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > stress_vm; ///< Multicomponent field variable storing pointwise strains (6 components)
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > body_force; ///< Multicomponent field variable storing pointwise displacements
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > energy; ///< Single component variable storing strain energy
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > energies; ///< Multicomponent field variable storing individual strain energies

	BC::BC *mybc;

	std::string ic_type, filename;

	IC::IC *ic;

	struct
	{
		bool        on       = false;
		Grid        grid = Grid::Node;
		int         interval = 1;
		int         max_iter = 20;
		int         verbose = 3;
		int         cgverbose = 3;
		amrex::Real tol_rel = 0.0;
		amrex::Real tol_abs = 1.0E-10;
		amrex::Real tstart = 0.0;
		amrex::Vector<amrex::Real> load_t;
		amrex::Vector<amrex::Real> load_disp;

		// Node only
		amrex::Vector<amrex::Real> AMREX_D_DECL(bc_xlo,bc_ylo,bc_zlo);
		amrex::Vector<amrex::Real> AMREX_D_DECL(bc_xhi,bc_yhi,bc_zhi);
		std::array<Operator::Elastic<model_type>::BC,AMREX_SPACEDIM> AMREX_D_DECL(bctype_xlo, bctype_ylo, bctype_zlo);
		std::array<Operator::Elastic<model_type>::BC,AMREX_SPACEDIM> AMREX_D_DECL(bctype_xhi, bctype_yhi, bctype_zhi);


		Operator::Elastic<model_type> *op;
		//BC::BC *bc;

	} elastic;
	//#endif

	BC::BC *elastic_bc;

	amrex::Vector<amrex::FabArray<amrex::BaseFab<model_type> > > model;
	model_type modeltype;
	Set::Scalar rho = 1.0;
	std::string 	input_material = "isotropic";

	const std::array<amrex::IntVect,AMREX_SPACEDIM> dx = {{AMREX_D_DECL(amrex::IntVect(AMREX_D_DECL(1,0,0)),
									    amrex::IntVect(AMREX_D_DECL(0,1,0)),
									    amrex::IntVect(AMREX_D_DECL(0,0,1)))}};

};
}
#endif
