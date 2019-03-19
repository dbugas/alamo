#include "ElasticDynamics.H"
#include "BC/Constant.H"
#include "Set/Set.H"
#include "Util/Util.H"
#include "IC/Random.H"
#include "IC/Trig.H"
#include "Model/Solid/LinearElastic/Isotropic.H"

namespace Integrator
{
ElasticDynamics::ElasticDynamics() : Integrator()
{
	//
	// READ INPUT PARAMETERS
	//
	{
		amrex::ParmParse pp("amr");
		pp.query("max_level",max_level);
	}

	{
		amrex::ParmParse pp_material("material");
		pp_material.query("model",input_material);
		if(input_material == "isotropic")
		{
			amrex::Real lambda = 3.0;
			amrex::Real mu = 4.0;
			amrex::ParmParse pp_material_isotropic("material.isotropic");
			pp_material_isotropic.query("lambda",lambda);
			pp_material_isotropic.query("mu",mu);
			pp_material_isotropic.query("rho",rho);
			if(lambda <=0)
			{
				Util::Warning(INFO, "Warning. Lambda must be positive. Resetting back to default value");
				lambda = 1.0;
			}
			if(mu <= 0)
			{
				Util::Warning(INFO, "Warning. Mu must be positive. Resetting back to default value");
				mu = 1.0;
			}
			modeltype = model_type(lambda,mu);
			//modeltype = model_type(1.0);
			//Model::Solid::Elastic::Degradable::Isotropic modeltype(lambda,mu);
			//modeltype(lambda,mu);
		}
		else if(input_material == "cubic")
			Util::Abort(INFO, "Not implemented yet");
	}

	// Elasticity
	{
		amrex::ParmParse pp("elastic");
		pp.query("on",elastic.on);

		if (elastic.on)
		{
			pp.query("interval",elastic.interval);
			pp.query("max_iter",elastic.max_iter);
			pp.query("verbose",elastic.verbose);
			pp.query("cgverbose",elastic.cgverbose);
			pp.query("tol_rel",elastic.tol_rel);
			pp.query("tol_abs",elastic.tol_abs);
			pp.query("tstart",elastic.tstart);

			std::string elastic_grid_str;
			pp.query("grid",elastic_grid_str);
			if (elastic_grid_str == "node") elastic.grid = Grid::Node;
			else if (elastic_grid_str == "cell") Util::Abort(INFO, "Cell elasitc model is not being used here.");


			if (elastic.grid == Grid::Node)
			{
				amrex::ParmParse pp_bc("elastic.bc");

				// Read in boundary types as strings, then convert to Operator::Elastic::BC types and store for future use.
				amrex::Vector<std::string> AMREX_D_DECL(bctype_xlo_str,bctype_ylo_str,bctype_zlo_str);
				amrex::Vector<std::string> AMREX_D_DECL(bctype_xhi_str,bctype_yhi_str,bctype_zhi_str);
				
				AMREX_D_TERM(pp_bc.queryarr("type_xlo",bctype_xlo_str);, pp_bc.queryarr("type_ylo",bctype_ylo_str);, pp_bc.queryarr("type_zlo",bctype_zlo_str););
				AMREX_D_TERM(pp_bc.queryarr("type_xhi",bctype_xhi_str);, pp_bc.queryarr("type_yhi",bctype_yhi_str);, pp_bc.queryarr("type_zhi",bctype_zhi_str););
				
				if ( AMREX_D_TERM(bctype_xlo_str.size() < AMREX_SPACEDIM, || bctype_ylo_str.size() < AMREX_SPACEDIM, || bctype_zlo_str.size() < AMREX_SPACEDIM)
			     	||
			     	AMREX_D_TERM(bctype_xhi_str.size() < AMREX_SPACEDIM, || bctype_yhi_str.size() < AMREX_SPACEDIM, || bctype_zhi_str.size() < AMREX_SPACEDIM))
					Util::Abort(INFO, "incorrect number of terms specified in bctype");

				std::map<std::string,Operator::Elastic<model_type>::BC> bc;
				bc["displacement"] = Operator::Elastic<model_type>::BC::Displacement;
				bc["disp"]         = Operator::Elastic<model_type>::BC::Displacement;
				bc["traction"]     = Operator::Elastic<model_type>::BC::Traction;
				bc["trac"]         = Operator::Elastic<model_type>::BC::Traction;
				bc["periodic"]     = Operator::Elastic<model_type>::BC::Periodic;
				
				AMREX_D_TERM(elastic.bctype_xlo = {AMREX_D_DECL(bc[bctype_xlo_str[0]], bc[bctype_xlo_str[1]], bc[bctype_xlo_str[2]])};,
				     elastic.bctype_ylo = {AMREX_D_DECL(bc[bctype_ylo_str[0]], bc[bctype_ylo_str[1]], bc[bctype_ylo_str[2]])};,
				     elastic.bctype_zlo = {AMREX_D_DECL(bc[bctype_zlo_str[0]], bc[bctype_zlo_str[1]], bc[bctype_zlo_str[2]])};);
				AMREX_D_TERM(elastic.bctype_xhi = {AMREX_D_DECL(bc[bctype_xhi_str[0]], bc[bctype_xhi_str[1]], bc[bctype_xhi_str[2]])};,
				     elastic.bctype_yhi = {AMREX_D_DECL(bc[bctype_yhi_str[0]], bc[bctype_yhi_str[1]], bc[bctype_yhi_str[2]])};,
				     elastic.bctype_zhi = {AMREX_D_DECL(bc[bctype_zhi_str[0]], bc[bctype_zhi_str[1]], bc[bctype_zhi_str[2]])};);


				AMREX_D_TERM(pp_bc.queryarr("xlo",elastic.bc_xlo);, pp_bc.queryarr("ylo",elastic.bc_ylo);, pp_bc.queryarr("zlo",elastic.bc_zlo););
				AMREX_D_TERM(pp_bc.queryarr("xhi",elastic.bc_xhi);, pp_bc.queryarr("yhi",elastic.bc_yhi);, pp_bc.queryarr("zhi",elastic.bc_zhi););


				RegisterNodalFab(displacement, AMREX_SPACEDIM,2,"disp");
				RegisterNodalFab(displacement_nm1, AMREX_SPACEDIM, 2, "disp_nm1");
				RegisterNodalFab(displacement_nm2, AMREX_SPACEDIM, 2, "disp_nm2");
				RegisterNodalFab(velocity, AMREX_SPACEDIM, 2, "vel");
				RegisterNodalFab(body_force, AMREX_SPACEDIM,2,"rhs");
				RegisterNodalFab(residual, AMREX_SPACEDIM,2,"resid");
				RegisterNodalFab(stress, AMREX_SPACEDIM*AMREX_SPACEDIM,2,"stress");
				RegisterNodalFab(divsigma, AMREX_SPACEDIM,2,"divSigma");
			}
		}
	}
}


void
ElasticDynamics::Advance (int lev, amrex::Real time, amrex::Real dt)
{

	/// TODO Make this optional
	if (lev != max_level) return;

	const amrex::Real* DX = geom[lev].CellSize();

	static amrex::IntVect AMREX_D_DECL(dx(AMREX_D_DECL(1,0,0)),
					   dy(AMREX_D_DECL(0,1,0)),
					   dz(AMREX_D_DECL(0,0,1)));

	for ( amrex::MFIter mfi(*displacement[lev],true); mfi.isValid(); ++mfi )
	{
		const amrex::Box& bx = mfi.validbox();
		FArrayBox &displacement_box     = (*displacement[lev])[mfi];
		FArrayBox &displacement_box_nm1     = (*displacement_nm1[lev])[mfi];
		FArrayBox &displacement_box_nm2     = (*displacement_nm2[lev])[mfi];
		FArrayBox &divsigma_box     = (*divsigma[lev])[mfi];
		FArrayBox &rhs_box     = (*body_force[lev])[mfi];
		FArrayBox &velocity_box = (*velocity[lev])[mfi];

		AMREX_D_TERM(for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++),
			     for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++),
			     for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++))
		{
			amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));
			for (int n = 0; n<AMREX_SPACEDIM; n++)
			{
				displacement_box(m,n) = 2.0*displacement_box_nm1(m,n) - displacement_box_nm2(m,n) + dt*dt*(divsigma_box(m,n) - rhs_box(m,n))/rho;
				velocity_box(m,n) = (displacement_box(m,n) - displacement_box_nm1(m,n))/dt;
			}
		}
	}
}

void
ElasticDynamics::Initialize (int lev)
{
	if (elastic.on)
	{
		displacement[lev].get()->setVal(0.0);
		displacement_nm1[lev].get()->setVal(0.0);
		displacement_nm2[lev].get()->setVal(0.0);
		velocity[lev].get()->setVal(0.0);
		strain[lev].get()->setVal(0.0); 
		stress[lev].get()->setVal(0.0); 
		divsigma[lev].get()->setVal(0.0);
		stress_vm[lev].get()->setVal(0.0);
		body_force[lev].get()->setVal(0.0);
		energy[lev].get()->setVal(0.0); 
		energies[lev].get()->setVal(0.0); 
	}
}


void
ElasticDynamics::TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
{
	const amrex::Real* DX      = geom[lev].CellSize();

	amrex::Vector<int>  itags;

	for (amrex::MFIter mfi(*displacement[lev],true); mfi.isValid(); ++mfi)
	{
		const amrex::Box&  bx  = mfi.validbox();
		amrex::TagBox&     tag  = tags[mfi];
		amrex::BaseFab<amrex::Real> &displacement_box = (*displacement[lev])[mfi];

		AMREX_D_TERM(	for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]-1; i++),
						for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]-1; j++),
						for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]-1; k++) )
		{
			for (int n = 0; n < AMREX_SPACEDIM; n++)
			{

				Set::Scalar cell_value = (1.0*(AMREX_D_TERM(2,+2,+4)))*(
											AMREX_D_TERM(displacement_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j,k))),
												+ displacement_box(amrex::IntVect(AMREX_D_DECL(i,j+1,k))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j+1,k))),
												+ displacement_box(amrex::IntVect(AMREX_D_DECL(i,j,k+1))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j,k+1)))
												+ displacement_box(amrex::IntVect(AMREX_D_DECL(i,j+1,k+1))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j+1,k+1)))));
				AMREX_D_TERM(
				Set::Scalar cell_value_left = (1.0*(AMREX_D_TERM(2,+2,+4)))*(
											AMREX_D_TERM(displacement_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i-1,j,k))),
												+ displacement_box(amrex::IntVect(AMREX_D_DECL(i,j+1,k))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i-1,j+1,k))),
												+ displacement_box(amrex::IntVect(AMREX_D_DECL(i,j,k+1))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i-1,j,k+1)))
												+ displacement_box(amrex::IntVect(AMREX_D_DECL(i,j+1,k+1))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i-1,j+1,k+1)))));

				Set::Scalar cell_value_right = (1.0*(AMREX_D_TERM(2,+2,+4)))*(
											AMREX_D_TERM(displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j,k))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+2,j,k))),
												+ displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j+1,k))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+2,j+1,k))),
												+ displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j,k+1))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+2,j,k+1)))
												+ displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j+1,k+1))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+2,j+1,k+1)))));
				,
				Set::Scalar cell_value_top = (1.0*(AMREX_D_TERM(2,+2,+4)))*(
											AMREX_D_TERM(displacement_box(amrex::IntVect(AMREX_D_DECL(i,j+1,k))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j+1,k))),
												+ displacement_box(amrex::IntVect(AMREX_D_DECL(i,j+2,k))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j+2,k))),
												+ displacement_box(amrex::IntVect(AMREX_D_DECL(i,j+1,k))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j+1,k+1)))
												+ displacement_box(amrex::IntVect(AMREX_D_DECL(i,j+2,k+1))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j+2,k+1)))));

				Set::Scalar cell_value_bottom = (1.0*(AMREX_D_TERM(2,+2,+4)))*(
											AMREX_D_TERM(displacement_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j,k))),
												+ displacement_box(amrex::IntVect(AMREX_D_DECL(i,j-1,k))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j-1,k))),
												+ displacement_box(amrex::IntVect(AMREX_D_DECL(i,j-1,k))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j,k+1)))
												+ displacement_box(amrex::IntVect(AMREX_D_DECL(i,j-1,k+1))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j-1,k+1)))));
				,
				Set::Scalar cell_value_back = (1.0*(AMREX_D_TERM(2,+2,+4)))*(
											AMREX_D_TERM(displacement_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j,k))),
												+ displacement_box(amrex::IntVect(AMREX_D_DECL(i,j+1,k))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j+1,k))),
												+ displacement_box(amrex::IntVect(AMREX_D_DECL(i,j,k-1))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j,k-1)))
												+ displacement_box(amrex::IntVect(AMREX_D_DECL(i,j+1,k-1))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j+1,k-1)))));

				Set::Scalar cell_value_front = (1.0*(AMREX_D_TERM(2,+2,+4)))*(
											AMREX_D_TERM(displacement_box(amrex::IntVect(AMREX_D_DECL(i,j,k+1))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j,k+1))),
												+ displacement_box(amrex::IntVect(AMREX_D_DECL(i,j+1,k+1))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j+1,k+1))),
												+ displacement_box(amrex::IntVect(AMREX_D_DECL(i,j,k+2))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j,k+2)))
												+ displacement_box(amrex::IntVect(AMREX_D_DECL(i,j+1,k+2))) + displacement_box(amrex::IntVect(AMREX_D_DECL(i+1,j+1,k+2)))));
				);

				AMREX_D_TERM(	Set::Scalar gradx = (cell_value_right - cell_value_left)/(2.0*DX[0]); ,
								Set::Scalar grady= (cell_value_top - cell_value_bottom)/(2.0*DX[1]); ,
								Set::Scalar gradz = (cell_value_front - cell_value_back)/(2.0*DX[2]);
				);
				if (DX[0]*sqrt(AMREX_D_TERM(gradx*gradx, + grady*grady, + gradz*gradz)) > 0.1) tag(amrex::IntVect(AMREX_D_DECL(i,j,k))) = amrex::TagBox::SET;
			}
		}
	}
}


void ElasticDynamics::TimeStepComplete(amrex::Real /*time*/, int /*iter*/)
{
	int nlevels = maxLevel() + 1;
	for (int ilev = 0; ilev<nlevels; ilev++)
	{
		amrex::MultiFab::Copy(*displacement_nm2[ilev], *displacement_nm1[ilev], 0, 0, AMREX_SPACEDIM, 2);
		amrex::MultiFab::Copy(*displacement_nm1[ilev], *displacement[ilev], 0, 0, AMREX_SPACEDIM, 2);
	}
}


void ElasticDynamics::TimeStepBegin(amrex::Real time, int iter)
{
	if (!elastic.on) return;
	if (iter%elastic.interval) return;
	if (time < elastic.tstart) return;

	if (elastic.grid == Grid::Node)
	{

		LPInfo info;
		info.setAgglomeration(true);
		info.setConsolidation(true);
		int max_mg_level = 0;
		info.setMaxCoarseningLevel(max_mg_level);

		int nlevels = maxLevel() + 1;
		amrex::Vector<amrex::BoxArray> 	ngrids;
		amrex::Vector<amrex::FabArray<amrex::BaseFab<model_type> > > model;

		ngrids.resize(nlevels);
		model.resize(nlevels);

		//IC::Trig ic(geom);
		//std::complex<int> I(0,1);
		//ic.Define(Util::Random(),AMREX_D_DECL(I,I,I));

		for (int ilev = 0; ilev < nlevels; ++ilev)
		{
			ngrids[ilev] = grids[ilev];
			ngrids[ilev].convert(amrex::IntVect::TheNodeVector());

			model[ilev].define(ngrids[ilev],dmap[ilev],1,2);

			Util::Message(INFO,"displacement size = " , displacement.size());
			Util::Message(INFO,"body_force size = " , body_force.size());
			
			displacement[ilev]->setVal(0.0);
			stress[ilev]->setVal(0.0);
			body_force[ilev]->setVal(0.000000001);

			model[ilev].setVal(modeltype);

			
			for (amrex::MFIter mfi(*body_force[ilev],true); mfi.isValid(); ++mfi)
			{
			 	const amrex::Box& box = mfi.validbox();

				amrex::BaseFab<amrex::Real> &rhsfab = (*(body_force[ilev]))[mfi];
				AMREX_D_TERM(for (int i = box.loVect()[0]; i<=box.hiVect()[0]+1; i++),
				 	     for (int j = box.loVect()[1]; j<=box.hiVect()[1]+1; j++),
				 	     for (int k = box.loVect()[2]; k<=box.hiVect()[2]+1; k++))
				{
					amrex::IntVect m(AMREX_D_DECL(i,j,k));

					for (int p = 0; p<AMREX_SPACEDIM; p++)
					{
						AMREX_D_TERM( if (i == geom[ilev].Domain().loVect()[0]) rhsfab(m,p)   = elastic.bc_xlo[p];,
							      if (j == geom[ilev].Domain().loVect()[1]) rhsfab(m,p)   = elastic.bc_ylo[p];,
							      if (k == geom[ilev].Domain().loVect()[2]) rhsfab(m,p)   = elastic.bc_zlo[p]; );
						AMREX_D_TERM( if (i == geom[ilev].Domain().hiVect()[0]+1) rhsfab(m,p) = elastic.bc_xhi[p];,
							      if (j == geom[ilev].Domain().hiVect()[1]+1) rhsfab(m,p) = elastic.bc_yhi[p];,
							      if (k == geom[ilev].Domain().hiVect()[2]+1) rhsfab(m,p) = elastic.bc_zhi[p]; );
					}

				}
			}
		}
		
		elastic.op = new Operator::Elastic<model_type>();
		elastic.op->define(geom,grids,dmap,info);

		elastic.op->SetBC({{AMREX_D_DECL(elastic.bctype_xlo,elastic.bctype_ylo,elastic.bctype_zlo)}},
				  {{AMREX_D_DECL(elastic.bctype_xhi,elastic.bctype_yhi,elastic.bctype_zhi)}});

		for (int ilev = 0; ilev < nlevels; ++ilev) elastic.op->SetModel(ilev,model[ilev]);

		amrex::MLMG solver(*elastic.op);
		solver.apply(GetVecOfPtrs(divsigma),GetVecOfPtrs(displacement_nm1));
	}
}

}
