/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include "pair_POL.h"
#include "atom.h"
#include "molecule.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "fix_bd_baoab.h"
#include "random_mars.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

PairPOL::PairPOL(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairPOL::~PairPOL()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(Jmin);
    memory->destroy(rad);
    memory->destroy(k_at);
    memory->destroy(cut);
    memory->destroy(fmag);
    memory->destroy(dotsign);
    memory->destroy(offset);
	memory->destroy(PolNeigh);
	memory->destroy(NpolNeigh);
	memory->destroy(pollist);
	memory->destroy(f_rand);
  }
}

/* ---------------------------------------------------------------------- */

void PairPOL::compute(int eflag, int vflag)
{
	int i,j,ii,jj,inum,jnum,itype,jtype, c_ind;
	double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair, dx, dy, dz, dnorm;
	double fvecx, fvecy, fvecz, fvecnorm, dxhat, dyhat, dzhat, delnorm, dx0, dx1, dy0, dy1, dz0, dz1;
	double rsq,r2inv,r,rinv,factor, dotprod;
	int *ilist,*jlist,*numneigh,**firstneigh;
	//std::cout << "In compute\n";
	evdwl = 0.0;
	if (eflag || vflag) ev_setup(eflag,vflag);
	else evflag = vflag_fdotr = 0;

	double **x = atom->x;
	double **f = atom->f;
	int *type = atom->type;
	int nlocal = atom->nlocal;
	double *special_lj = force->special_lj;
	int newton_pair = force->newton_pair;
	double * dx_mol;

	// Set up arrays for lists
	inum = list->inum;
	ilist = list->ilist;
	numneigh = list->numneigh;
	firstneigh = list->firstneigh;

	// Create array of polymerase indices
	int n_pol=0;
	for (ii = 0; ii < inum; ii++){
		i = ilist[ii];
		if (type[i]!=chain_type){
			pollist[n_pol] = i;
			n_pol++;
			if(n_pol>Npol){
				//std::cout << "Error: n_pol exceeds Npol: " << n_pol << "\n";
				exit(1);
			}
		}
	}
	// Initialize NpolNeigh (array of number of neighbors for each polymerase)
	//std::cout << "There are "<< n_pol << "polymerases, Npol: " << Npol << " \n";
	for(int ind=0; ind<Npol; ind++){
		NpolNeigh[ind]=0;
	}
	/*
	for (ii = 0; ii < inum; ii++) {
		i = ilist[ii];
		if (type[i]!=chain_type){
			pollist[n_pol] = i;
			n_pol++;
		}
	}
	*/
	//Fill NpolNeigh (Array of number of neighbors for each polymerase) and PolNeigh (Array of polymerase neighbor list arrays)
	for (ii = 0; ii < inum; ii++) {
		i = ilist[ii];
		jlist = firstneigh[i];
	    jnum = numneigh[i];
	    //std::cout << "atom " << atom->tag[i] << " has " << jnum << " neighbors : \n\t"; 
		for(jj=0; jj<jnum; jj++){
			j = jlist[jj];
			factor = special_lj[sbmask(j)];
			j &= NEIGHMASK;
			for(int kk=0; kk<Npol; kk++){
				if((pollist[kk]==i) || (pollist[kk]==j)){
					if(type[i]==chain_type){
						PolNeigh[kk][NpolNeigh[kk]]=i;
						NpolNeigh[kk]++;
					}else{
						PolNeigh[kk][NpolNeigh[kk]]=j;
						NpolNeigh[kk]++;
					}
				}
			}
			//std:: cout << atom->tag[j] << " ";
		}
		//std::cout << "\n";
	}
	/*
	std::cout << "Neighbor list for "<< Npol << " Polymerases \n";
	// Check neighbor list for pols
	for (ii = 0; ii < Npol; ii++) {
		i = pollist[ii];
	    std::cout << "\natom " << atom->tag[i] << " has " << NpolNeigh[ii] << " neighbors : \n\t"; 
		for(jj =0; jj<NpolNeigh[ii]; jj++){
			j = PolNeigh[ii][jj];
			std:: cout << atom->tag[j] << " ";
		}
	}
	std::cout << "\n";
	*/
	//std::cout << "About to make polymerase updates \n";
	// Loop through all polymerases
	for(ii = 0; ii<Npol; ii++){
		//std::cout << "ii: " << ii <<", Polymerase id: " << atom->tag[pollist[ii]] << "\n";  
		if(NpolNeigh[ii]>=1){
			i = pollist[ii];

			xtmp = x[i][0];
			ytmp = x[i][1];
			ztmp = x[i][2];
			itype = type[i];
			
			j=PolNeigh[ii][0];
			delx = xtmp - x[j][0];
			dely = ytmp - x[j][1];
			delz = ztmp - x[j][2];
			double r_min = delx*delx + dely*dely + delz*delz;
			int j_min = j;	

			// Find chain index closest to polymerase
			for(jj = 1; jj < NpolNeigh[ii]; jj++){
				j=PolNeigh[ii][jj];

				delx = xtmp - x[j][0];
				dely = ytmp - x[j][1];
				delz = ztmp - x[j][2];
				rsq = delx*delx + dely*dely + delz*delz;
				if(rsq<r_min){
					r_min = rsq;
					j_min = j;
				}
			}
			Jmin[ii] = j_min; // Log index of closest j to polymerase ii
			//std::cout << "ID of nearest " <<  atom->tag[j_min] << ", boundlow: " << boundlow << "\n";
			
			if (atom->tag[j_min] == boundlow){
				// Reassign Pol
				/*
				std::cout << "Reassigning " << atom->tag[i] << "\n\n";
				std::cout << " boundlow: " << boundlow << ", boundhigh: " << boundhigh << "\n";
				std::cout << "Before: x[" << atom->tag[i] << "]: " << x[i][0] << ", " << x[i][1] << ", " << x[i][2] << "\n";
				*/
				int bound_local = atom->map(boundhigh);
				int bound_local2 = atom->map(boundhigh-1);
				double x_new = 0.5*(x[bound_local][0] + x[bound_local2][0]);
				double y_new = 0.5*(x[bound_local][1] + x[bound_local2][1]);
				double z_new = 0.5*(x[bound_local][2] + x[bound_local2][2]);
				double dist2 = sigma2 +1;
				for(int pol_ind=0; pol_ind<Npol; pol_ind++){
					if(pol_ind!=ii){
						int ind_p = pollist[pol_ind];
						double dist0 = (x[ind_p][0] - x_new) * (x[ind_p][0] - x_new) +  (x[ind_p][1] - y_new) * (x[ind_p][1] - y_new) + (x[ind_p][2] - z_new) * (x[ind_p][2] - z_new);
						//std::cout << "pol: " << atom->tag[ind_p] <<", r: " << x[ind_p][0] << " " << x[ind_p][1] << " " << x[ind_p][2] << "\n";
						//std::cout << "dist0: " << dist0 << "\n";
						if(dist0<dist2){
							dist2 = dist0;
						}
					}
				}
				if (dist2>sigma2){
					x[i][0] = x_new;
					/*
					std::cout << "dist2: " << dist2 << ", sigma2:" << sigma2 << "\n";
					std::cout << "Bound_local: " << atom->tag[bound_local] << "\n";
					std::cout << "Bound_local - 1: " << atom->tag[bound_local2] << "\n";
					std::cout << 0.5*(x[bound_local][0] + x[bound_local2][0]) << "\n";
					std::cout << 0.5*(x[bound_local][1] + x[bound_local2][1]) << "\n";
					std::cout << 0.5*(x[bound_local][2] + x[bound_local2][2]) << "\n";
					*/
					x[i][1] = y_new;
					x[i][2] = z_new;
					f[i][0]=0.;
					f[i][1]=0.;
					f[i][2]=0.;
					Jmin[ii] = bound_local;
				}
				/*
				else{
					std::cout << "Bound atom is currently occupied\n";
					std::cout << "dist2: " << dist2 << ", sigma2:" << sigma2 << "\n";
				}*/
				
				//std::cout << "x[" << atom->tag[i] << "]: " << x[i][0] << ", " << x[i][1] << ", " << x[i][2] << "\n";
				//exit(1);
			}//else{

			// Add Force Calculations

			j = j_min;

			delx = xtmp - x[j][0];
			dely = ytmp - x[j][1];
			delz = ztmp - x[j][2];
			rsq = delx*delx + dely*dely + delz*delz;
		    jtype = type[j];

		  	//std::cout << "B\n";
		  	if (rsq < cutsq[itype][jtype]) {

			    //std::cout << "Calculating Pol Potential and Force\n";
			    //r2inv = 1.0/rsq;
			    //r = sqrt(rsq);
			    //rinv = 1.0/r;
			    
			    /*  
				//std::cout << "C\n";
				std::cout << "atom " << i << " belongs to molecule: " 
					<< atom->molecule[i] << " and has type: " << type[i] << std::endl;
				std::cout << "atom " << j << " belongs to molecule: " 
					<< atom->molecule[j] << " and has type: " << type[j] << std::endl;
			    tagint* bond_atom_i = atom->bond_atom[i];
			    int nbond_i = atom->num_bond[i];
				std::cout << "atom " << i << " has " << nbond_i << " bonds: ";
			    for(int indval=0; indval<nbond_i; indval++){
			    	std::cout<< bond_atom_i[indval] << " ";
			    }
			    std::cout << std::endl;

			    bond_atom_i = atom->special[i];
			    nbond_i = atom->nspecial[i][0];
				std::cout << "atom " << i << " has " << nbond_i << " bonds: ";
			    for(int indval=0; indval<nbond_i; indval++){
			    	std::cout<< bond_atom_i[indval] << " " << " mapping: " << atom->map(bond_atom_i[indval]) << " ";
			    }
			    std::cout << std::endl;
			    
			    tagint* bond_atom_j = atom->bond_atom[j];
			    int nbond_j = atom->num_bond[j];
				std::cout << "atom " << j << " has " << nbond_j << " bonds: ";
			    for(int indval=0; indval<nbond_j; indval++){
			    	std::cout<< bond_atom_j[indval] << " ";
			    }
			    std::cout << std::endl;
			    bond_atom_j = atom->special[j];
			    nbond_j = atom->nspecial[j][0];
				std::cout << "atom " << j << " with id: " << atom->tag[j] << " check: " << atom->map(atom->tag[j]) 
					<< " has " << nbond_j << " bonds: ";
			    for(int indval=0; indval<nbond_j; indval++){
			    	std::cout<< bond_atom_j[indval] << " "<< " mapping: " << atom->map(bond_atom_j[indval]) << ", ";
			    }
			    std::cout << std::endl;
			    
			    //std::cout << ", dx_ind: " << dx_ind[itype][jtype] << "\n";
			    */
			    // have force point from high ind to low ind
			    int ind_mol0 = 0;
			    int c_ind = 0;
			    if(dotsign[itype][jtype]>0){
			    	ind_mol0=i; //Index of pol atom
			    	c_ind = j; //Index of chain monomer
			    }else{
			    	ind_mol0=j; //Index of pol atom
			    	c_ind = i; // index of chain monomer
			    }
			    //int c0 = atom->map(atom->special[c_ind][0]);
			    //int c1 = atom->map(atom->special[c_ind][1]);
			    //std::cout << " Attempt to get neighboring atoms along chain \n";
			    int c0 = atom->map(atom->tag[c_ind]-1);
			    int c1 = atom->map(atom->tag[c_ind]+1);
			    if(atom->tag[j_min] == boundhigh){
			    	c1=c0;
			    }
			    if(c0<0 || c1<0){
			    	std::cout << "Error: c0 or c1 not proper \n";
			    	std:: cout << "  c1: " << c1 << ", c0: " << c0 << "\n";
			    	std:: cout << "  c_ind: " << c_ind << ", id: " << atom->tag[c_ind] << ", nlocal: " << atom->nlocal << "\n";
			    	return;
			    }
			    //std::cout << "Got em! Calculating displacement vectors ... \n";
			    /*
			    if((atom->nspecial[c_ind][0])<1){
			      std::cout << "Calculate dr0 and dr1 \n";
			      int ind = atom->tag[c_ind];
			      ind ++;
			      std::cout << atom->tag[c_ind] << " has " << atom->nspecial[c_ind][0] << " bonds\n";
			      std::cout << ind << " maps to " << atom->map(ind) << "\n";
			      ind -=2;
			      std::cout <<  ind << " maps to " << atom->map(ind) << "\n";
			      std::cout << "c0: " << c0 << " has id: " << atom->tag[c0] <<"\n c1:" << c1  << " has id: " << atom->tag[c0] << "\n";
			    }
			    */
			    // Create vectors pointing from the POL atom to the chain neighbors of the interacting monomer
			    dx0 = x[c0][0] - x[ind_mol0][0];
			    //std::cout << "Calculated dx0 \n" << x[c0][0];
			    dy0 = x[c0][1] - x[ind_mol0][1];
			    dz0 = x[c0][2] - x[ind_mol0][2];
			    dx1 = x[c1][0] - x[ind_mol0][0];
			    //std::cout << "Calculated dx1 \n" << x[c1][0];
			    dy1 = x[c1][1] - x[ind_mol0][1];
			    dz1 = x[c1][2] - x[ind_mol0][2];

			    double x10=0.;
			    double y10=0.;
			    double z10=0.;
			    
			    //std::cout << "Calculate r10\n";
			    if((dx0*dx0+dy0*dy0+dz0*dz0)<(dx1*dx1+dy1*dy1+dz1*dz1)){
			      // Pol atom is closer to c0 (c0 is r1)
			      x10= x[c0][0] - x[c_ind][0];
			      y10= x[c0][1] - x[c_ind][1];
			      z10= x[c0][2] - x[c_ind][2];
			      c1=c0;
			    }else{
			      // Pol atom is closer to c1
			      x10= x[c1][0] - x[c_ind][0];
			      y10= x[c1][1] - x[c_ind][1];
			      z10= x[c1][2] - x[c_ind][2];
			    }

			    //std::cout << "Calculatd r10\n";
			    dx = x10;
			    dy = y10;
			    dz = z10;

			    if( (atom->tag[c1]) > (atom->tag[c_ind]) ){
			      // c1 has a higher index that c_ind, dx should point from c1 to c_ind
			      dx = -dx;
			      dy = -dy;
			      dz = -dz; 
			    } // else, dx = x10 etc., which we have already done
			    /*
			    std::cout << " Pol Atom: " << atom->tag[ind_mol0] << " is overlapping with chain atom " << atom->tag[c_ind] << 
			        ", and is close to chain atom " << atom->tag[c1] << "\n The vector from " << atom->tag[c_ind] << 
			        " to " << atom->tag[c1] << " is " << x10 << ", " << y10 << ", " << z10 << 
			        "\n The pulling force on pol atom points along " << dx << ", " << dy << ", " << dz <<"\n";
			    */
			    //std::cout << "Calculate r_hat\n";
			    dnorm = sqrt(dx*dx+dy*dy+dz*dz);
			    dxhat = dx/dnorm;
			    dyhat = dy/dnorm;
			    dzhat = dz/dnorm;

			    //std::cout << "Calculate l\n";
			    
			    double s = dotsign[itype][jtype]*(delx*x10 + dely*y10 + delz*z10)/(x10*x10+y10*y10+z10*z10);
			    if(s<0.){
			    	s=0.;
			    }
			    double lx = dotsign[itype][jtype]*delx - s*x10;
			    double ly = dotsign[itype][jtype]*dely - s*y10;
			    double lz = dotsign[itype][jtype]*delz - s*z10;
			    //double lsq = lx*lx + ly*ly +lz*lz;
			    //double l  = sqrt(lsq);
			    /*
			    std::cout << " The force of attraction for the pol atom points along " << -lx << ", " << -ly << "," << -lz
			      << "\n The dot product of this and the pulling force is " << lx*dx + ly*dy + lz*dz << "\n";
			    */
			    // Magnitude of pulling force
			    fpair = factor*fmag[itype][jtype];

			    // Generate random unit vector
			    double rx = random->gaussian();
        		double ry = random->gaussian();
        		double rz = random->gaussian();

        		double randnorm = sqrt(rx*rx + ry*ry + rz*rz);
        		rx /= randnorm;
        		ry /= randnorm;
        		rz /= randnorm;

			    //std::cout << "Calculate and assign forces\n";
			    //std::cout << " lsq: " << lsq << ", k_at: " << k_at[itype][jtype] << " factor: " << factor << "\n";
			    // Assign attractive interaction between line segment and point and pulling force
			    f[ind_mol0][0] += - factor * k_at[itype][jtype] * lx + fpair*dxhat + f_rand[itype][jtype]*rx;
			    f[ind_mol0][1] += - factor * k_at[itype][jtype] * ly + fpair*dyhat + f_rand[itype][jtype]*ry;
			    f[ind_mol0][2] += - factor * k_at[itype][jtype] * lz + fpair*dzhat + f_rand[itype][jtype]*rz;
			    /*
			    std::cout << "Local polymerase ind: " << ind_mol0 << " has id " << atom->tag[ind_mol0] << " and has attractive force: " << - factor * k_at[itype][jtype] * lsq * l * lx
						<< ", " << - factor * k_at[itype][jtype] * lsq * l * ly << ", " << - factor * k_at[itype][jtype] * lsq * l * lz << " and has total active force: "
			    		<< - factor * k_at[itype][jtype] * lsq * l * lx + fpair*dxhat << ", " << - factor * k_at[itype][jtype] * lsq * l * ly + fpair*dyhat <<
			    		", " << - factor * k_at[itype][jtype] * lsq * l * lz + fpair*dzhat << "\n\t It has total force: " << f[ind_mol0][0] << ", " << f[ind_mol0][1]
			    		<< ", " << f[ind_mol0][2] << "\n";
			    */
			    double sp = 0.5;
			    f[c1][0] += factor * k_at[itype][jtype] * s * lx - sp*fpair*dxhat + sp*f_rand[itype][jtype]*rx;
			    f[c1][1] += factor * k_at[itype][jtype] * s * ly - sp*fpair*dyhat + sp*f_rand[itype][jtype]*ry;
			    f[c1][2] += factor * k_at[itype][jtype] * s * lz - sp*fpair*dzhat + sp*f_rand[itype][jtype]*rz;
				/*
				std::cout << "Local chain ind: " << c1 << " has id " << atom->tag[c1] << " and has attractive force: " << factor * k_at[itype][jtype] * lsq * l * s * lx
						<< ", " << factor * k_at[itype][jtype] * lsq * l * s * ly << ", " << factor * k_at[itype][jtype] * lsq * l * s * lz << " and has total active force: "
			    		<< factor * k_at[itype][jtype] * lsq * l * s * lx - 0.5*fpair*dxhat << ", " << factor * k_at[itype][jtype] * lsq * l * s * ly - 0.5*fpair*dyhat <<
			    		", " << factor * k_at[itype][jtype] * lsq * l * s * lz - 0.5*fpair*dzhat << "\n\t It has total force: " << f[c1][0] << ", " << f[c1][1]
			    		<< ", " << f[c1][2] << "\n";
			    */
			    f[c_ind][0] += factor * k_at[itype][jtype] * (1-s) * lx - (1-sp)*fpair*dxhat + sp*f_rand[itype][jtype]*rx;
			    f[c_ind][1] += factor * k_at[itype][jtype] * (1-s) * ly - (1-sp)*fpair*dyhat + sp*f_rand[itype][jtype]*ry;
			    f[c_ind][2] += factor * k_at[itype][jtype] * (1-s) * lz - (1-sp)*fpair*dzhat + sp*f_rand[itype][jtype]*rz;
				/*
				std::cout << "Local chain ind2: " << c_ind << " has id " << atom->tag[c_ind] << " and has attractive force: " << factor * k_at[itype][jtype] * lsq * l * (1-s) * lx
						<< ", " << factor * k_at[itype][jtype] * lsq * l * (1-s) * ly << ", " << factor * k_at[itype][jtype] * lsq * l * (1-s) * lz << " and has total active force: " 
			    		<< factor * k_at[itype][jtype] * lsq * l * (1-s) * lx - 0.5*fpair*dxhat << ", " << factor * k_at[itype][jtype] * lsq * l * (1-s) * ly - 0.5*fpair*dyhat <<
			    		", " << factor * k_at[itype][jtype] * lsq * l * (1-s) * lz - 0.5*fpair*dzhat << "\n\t It has total force: " << f[c_ind][0] << ", " << f[c_ind][1]
			    		<< ", " << f[c_ind][2] << "\n";
			    */
			    if (eflag) {
			      // dr = r_A0 dot r_10 is the distance along the bond of pol from c_ind
			      double dr = dotsign[itype][jtype] * ( delx * x10 + dely * y10 + delz * z10);
			      if( atom->tag[c_ind] > atom->tag[c1]){
			        // Pulling towards c1 not c_ind, need to subtract from length of bond (|r_10|)
			        dr = sqrt( x10 * x10 + y10 * y10 + z10 * z10) - dr;
			      }
			      // U = Fpull*dr + 0.5kl^2
			      evdwl = fmag[itype][jtype] * dr + 0.5 * k_at[itype][jtype] * (lx * lx + ly * ly + lz * lz) - offset[itype][jtype];
			      evdwl *= factor;
			    }

			    if (evflag) ev_tally(i,j,nlocal,newton_pair,
			                         evdwl,0.0,fpair,delx,dely,delz);
			      
			    
			    //std::cout << "Calculated Pol Potential and Force\n";
			}	 
		}
	}
	/*
	// loop over neighbors of my atoms
	for (ii = 0; ii < inum; ii++) {
		i = ilist[ii];
		*/
		/*
		Molecule *onemol = atom->molecules[atom->molecule[i]-1];
		if(onemol->natoms <10){
			std::cout << "atom " << i << " belongs to molecule: " << atom->molecule[i] << " and has type: " << type[i] << std::endl;
		    std::cout << "r: (" << x[i][0] << "," << x[i][1] << "," << x[i][2] << ")\n";
			std::cout << "molecule " << atom->molecule[i]+1 << " contains, " << onemol->natoms << " atoms of " << onemol->ntypes << " types with "<< onemol->nbonds << " bonds \n";    
		    std::cout << " leading order orientation vector is: " << onemol->inertia[0] << "," << onemol->inertia[1] << "," << onemol->inertia[2] << "\n";
		    std::cout << " leading order orientation vector is: " << onemol->ex[0] << "," << onemol->ex[1] << "," << onemol->ex[2] << "\n";
		    std::cout << " second leading order orientation vector is: " << onemol->ey[0] << "," << onemol->ey[1] << "," << onemol->ey[2] << "\n";
		    std::cout << " third leading order orientation vector is: " << onemol->ez[0] << "," << onemol->ez[1] << "," << onemol->ez[2] << "\n";
		    for(int indval=0; indval<onemol->natoms; indval++){
		    	std::cout << "type of atom " << indval << " in molecule is " <<onemol->type[indval] << ", r = " << onemol->x[indval][0] << ", " << onemol->x[indval][1] << ", "<< onemol->x[indval][2] << "\n";
		    	std::cout << " displacement vector: " << onemol->dx[indval][0] << " " << onemol->dx[indval][1] << " " << onemol->dx[indval][2] << std::endl;
		    }
		}
		
	    jlist = firstneigh[i];
	    jnum = numneigh[i];
	    std::cout << "atom " << atom->tag[i] << " has " << jnum << " neighbors : \n\t"; 
		for(jj=0;jj<jnum;jj++){
			j = jlist[jj];
			factor = special_lj[sbmask(j)];
			j &= NEIGHMASK;
			std:: cout << atom->tag[j] << " ";
		}
		std::cout <<"\n";
		xtmp = x[i][0];
		ytmp = x[i][1];
		ztmp = x[i][2];
		itype = type[i];
		if(itype != chain_type){
		    jlist = firstneigh[i];
		    jnum = numneigh[i];
		    //std::cout << "A\n";

		    jj=0;
			j = jlist[jj];
			factor = special_lj[sbmask(j)];
			j &= NEIGHMASK;

			delx = xtmp - x[j][0];
			dely = ytmp - x[j][1];
			delz = ztmp - x[j][2];
			rsq = delx*delx + dely*dely + delz*delz;
			jtype = type[j];
		    int j_min = jlist[0];
		    double r_min = rsq; // index of j with minimum 
		    for (jj = 1; jj < jnum; jj++) {
		      j = jlist[jj];
		      factor = special_lj[sbmask(j)];
		      j &= NEIGHMASK;

		      delx = xtmp - x[j][0];
		      dely = ytmp - x[j][1];
		      delz = ztmp - x[j][2];
		      rsq = delx*delx + dely*dely + delz*delz;
		      if(rsq<r_min){
		      	r_min=rsq;
		      	j_min=jj;
		      }
		    }
		    jj = j_min;
		    j = jlist[jj];
			factor = special_lj[sbmask(j)];
			j &= NEIGHMASK;

			   
	    }  
	}
	*/
	/*
	std::cout << "Attempt to free memory associated with polymerase neighbor lists \n";
	memory->destroy(pollist);
	std::cout << "Successfully freed pollist \n";
	memory->destroy(PolNeigh);
	std::cout << "Successfully freed PolNeigh \n";
	memory->destroy(NpolNeigh);
	std::cout << "Successfully freed NpolNeigh \n";
	exit(1);
	*/
	//std::cout << "Processed pair_POL\n";
	if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairPOL::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  int N = atom->nmax;
  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(pollist, n,"pair:pollist");
  memory->create(NpolNeigh, n, "pair:NpolNeigh");
  memory->create(Jmin, n, "pair:Jmin");
  memory->create(PolNeigh, n, 40, "pair:PolNeigh");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(rad,n+1,"pair:rad");
  memory->create(k_at,n+1,n+1,"pair:k_at");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(fmag,n+1,n+1,"pair:fmag");
  //memory->create(chain_type,1,1,"pair:chain_type");
  //memory->create(chain_type,n+1,n+1,"pair:chain_type");
  memory->create(dotsign,n+1,n+1,"pair:dotsign");
  memory->create(offset,n+1,n+1,"pair:offset");
  memory->create(f_rand, n+1, n+1, "pair:f_rand");

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairPOL::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = utils::numeric(FLERR,arg[0],false,lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairPOL::coeff(int narg, char **arg)
{
  if (narg < 11 || narg > 12)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();
  
  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);
  /*
  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],atom->ntypes,ilo,ihi, error);
  utils::bounds(FLERR,arg[1],atom->ntypes,jlo,jhi, error);
  */
  double f_one = utils::numeric(FLERR,arg[2],false,lmp);
  chain_type = (int) utils::numeric(FLERR,arg[3],false,lmp);
  //int ctype_one = force->numeric(FLERR,arg[3]);
  double k_one = utils::numeric(FLERR,arg[4],false,lmp);
  Npol = (int) utils::numeric(FLERR, arg[5],false,lmp);
  boundlow = (int) utils::numeric(FLERR, arg[6],false,lmp);
  boundhigh = (int) utils::numeric(FLERR, arg[7],false,lmp);
  sigma2  = utils::numeric(FLERR, arg[8],false,lmp);
  double fr_one = utils::numeric(FLERR,arg[9],false,lmp);
  seed = (int) utils::numeric(FLERR,arg[10],false,lmp); //seed for random number generator. integer

  random = new RanMars(lmp,seed + comm->me);
  std::cout << "sigma2: " << sigma2 << "\n";
  sigma2 = sigma2 * sigma2;
  std::cout << "sigma2: " << sigma2 << "\n";
  std::cout << "boundlow: "<< boundlow <<"\n";
  std::cout << "boundhigh:" << boundhigh << "\n";
  //exit(0);
  //boundlow +=1;
  boundhigh -=1;
  std::cout << "Npol: " << Npol << "\n";
  memory->grow(pollist, Npol,"pair:pollist");
  memory->grow(NpolNeigh, Npol, "pair:NpolNeigh");
  memory->grow(PolNeigh, Npol, 40, "pair:PolNeigh");

  double cut_one = cut_global;
  if (narg == 8) cut_one = utils::numeric(FLERR,arg[6],false,lmp);
  //std::cout << "k_one: " << k_one << "\n";
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      fmag[i][j] = f_one;
      f_rand[i][j] = fr_one;
      //chain_type[i][j] = ctype_one;
      k_at[i][j] = k_one;
      //std::cout << " k_at[" << i <<"][" << j << "]: " << k_at[i][j] <<"\n";
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }
  //std::cout << "extracted pair coefficients\n";
  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairPOL::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    fmag[i][j] = mix_energy(fmag[i][i],fmag[j][j],1.0,1.0);
    f_rand[i][j] = mix_energy(f_rand[i][i],f_rand[j][j],1.0,1.0);
    k_at[i][j] = mix_energy(k_at[i][i],k_at[j][j],1.0,1.0);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }
  //std::cout << "In init_one, make dotsign \n";
  // Set dotsign
  if(i==chain_type){
  	dotsign[i][j]=-1;
  }else{
  	dotsign[i][j]=1;
  }


  if (offset_flag && (cut[i][j] > 0.0)) {
    offset[i][j] = 0.0;
  } else offset[i][j] = 0.0;


  fmag[j][i] = fmag[i][j];
  f_rand[j][i] = f_rand[i][j];
  k_at[j][i] = k_at[i][j];
  offset[j][i] = offset[i][j];
  dotsign[j][i] = - dotsign[i][j];
  //std::cout << "completed init_one \n";
  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairPOL::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&fmag[i][j],sizeof(double),1,fp);
        fwrite(&chain_type,sizeof(int),1,fp);
        fwrite(&k_at[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairPOL::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&fmag[i][j],sizeof(double),1,fp);
  	      fread(&chain_type,sizeof(int),1,fp);
          fread(&k_at[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&fmag[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&chain_type,1,MPI_INT,0,world);
        MPI_Bcast(&k_at[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairPOL::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairPOL::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairPOL::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g\n",i,fmag[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairPOL::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %d %d %g %g\n",i,j,fmag[i][j],dotsign[i][j],k_at[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairPOL::single(int i, int j, int itype, int jtype, double rsq,
                          double factor_coul, double factor_lj,
                          double &fforce)
{
	double r,rinv,phi;
	double fx, fy, fz, fpair;
	double **x = atom->x;
	double delx, dely, delz, dx0, dy0, dz0, dx1, dy1, dz1, dx, dy, dz;
	int ind_mol0 = 0;
	int c_ind = 0;
	int c0, c1;
	double x10,y10, z10;
	//std::cout << "In single\n";
	if(dotsign[itype][jtype]>0){
		ind_mol0=i; //Index of pol atom
		c_ind = j; //Index of chain monomer
	}else{
		ind_mol0=j; //Index of pol atom
		c_ind = i; // index of chain monomer
	}
	//std::cout << "i,j: " << atom->tag[i] << ", " << atom->tag[j] << "\n";
	int ii=0;
	for(ii=0; ii<Npol; ii++){
		//std::cout << "pollist[" << ii << "]: " << atom->tag[pollist[ii]] << " Npol: " << Npol << "\n";
		if (ind_mol0==pollist[ii]){
			// Then check if j is the nearest index or one off from the nearest index
			int j_id = atom->tag[j];
			int j_min = Jmin[ii];
			int jmin_id = atom->tag[j_min];
			//std::cout << "indmol0: " << atom->tag[ind_mol0] << ", jmin: " << atom->tag[j_min] << "\n"; 
			delx = x[ind_mol0][0] - x[j_min][0];
			dely = x[ind_mol0][1] - x[j_min][1];
			delz = x[ind_mol0][2] - x[j_min][2];
			int c0 = atom->map(atom->tag[j_min]-1);
		    int c1 = atom->map(atom->tag[j_min]+1);
			
			dx0 = x[c0][0] - x[ind_mol0][0];
		    //std::cout << "Calculated dx0 \n" << x[c0][0];
		    dy0 = x[c0][1] - x[ind_mol0][1];
		    dz0 = x[c0][2] - x[ind_mol0][2];
		    dx1 = x[c1][0] - x[ind_mol0][0];
		    //std::cout << "Calculated dx1 \n" << x[c1][0];
		    dy1 = x[c1][1] - x[ind_mol0][1];
		    dz1 = x[c1][2] - x[ind_mol0][2];


		    if((dx0*dx0+dy0*dy0+dz0*dz0)<(dx1*dx1+dy1*dy1+dz1*dz1)){
				// Pol atom is closer to c0 (c0 is r1)
				x10= x[c0][0] - x[j_min][0];
				y10= x[c0][1] - x[j_min][1];
				z10= x[c0][2] - x[j_min][2];
				c1=c0;
		    }else{
				// Pol atom is closer to c1
				x10= x[c1][0] - x[j_min][0];
				y10= x[c1][1] - x[j_min][1];
				z10= x[c1][2] - x[j_min][2];
		    }
		    // Only calculate if i or j is c_ind or c1
		    if(c_ind==j_min || c_ind == c1){

				c_ind = j_min;
			    //std::cout << "Calculatd r10\n";
			    dx = x10;
			    dy = y10;
			    dz = z10;

			    if( (atom->tag[c1]) > (atom->tag[c_ind]) ){
					// c1 has a higher index that c_ind, dx should point from c1 to c_ind
					dx = -dx;
					dy = -dy;
					dz = -dz; 
			    } 

				double dnorm = sqrt(dx*dx+dy*dy+dz*dz);
				double dxhat = dx/dnorm;
				double dyhat = dy/dnorm;
				double dzhat = dz/dnorm;

				double s = (delx*x10 + dely*y10 + delz*z10)/(x10*x10+y10*y10+z10*z10);
				if(s<0.){
					s=0.;
				}
				double lx = delx - s*x10;
				double ly = dely - s*y10;
				double lz = delz - s*z10;

				// Magnitude of pulling force
				fpair = factor_lj * fmag[itype][jtype];
				double sp = 0.5;
				// Assign attractive interaction between line segment and point and pulling force
				if(i==c1){
					fx = factor_lj * k_at[itype][jtype] * s * lx - sp*fpair*dxhat;
				    fy = factor_lj * k_at[itype][jtype] * s * ly - sp*fpair*dyhat;
				    fz = factor_lj * k_at[itype][jtype] * s * lz - sp*fpair*dzhat;
				}else if(i==c_ind){
					fx = factor_lj * k_at[itype][jtype] * (1-s) * lx - (1-sp)*fpair*dxhat;
				    fy = factor_lj * k_at[itype][jtype] * (1-s) * ly - (1-sp)*fpair*dyhat;
				    fz = factor_lj * k_at[itype][jtype] * (1-s) * lz - (1-sp)*fpair*dzhat;
				}else if(j==c1){
					fx = -factor_lj * k_at[itype][jtype] * s * lx + sp*fpair*dxhat;
				    fy = -factor_lj * k_at[itype][jtype] * s * ly + sp*fpair*dyhat;
				    fz = -factor_lj * k_at[itype][jtype] * s * lz + sp*fpair*dzhat;
				}else{
					// j==c_ind
					//std::cout << "(j,c_ind): " << j << ", " << c_ind << "\n";
					fx = -factor_lj * k_at[itype][jtype] * (1-s) * lx + (1-sp)*fpair*dxhat;
				    fy = -factor_lj * k_at[itype][jtype] * (1-s) * ly + (1-sp)*fpair*dyhat;
				    fz = -factor_lj * k_at[itype][jtype] * (1-s) * lz + (1-sp)*fpair*dzhat;
				}
				//std::cout << "s: " << s << ", l: " << lx << ", " << ly << ", " << lz << "\n";
				//std::cout << "Pair Forces i j fx fy fz\n";
				//std::cout << atom->tag[i] << "\t" << atom->tag[j] << "\t" << fx << "\t" << fy  << "\t" << fz << "\n";
				/*
				std::cout << atom->tag[i] << ", " << atom->tag[j] << ", f = " << fx << ", " << fy << ", " << fz 
					<< " c_ind: " << atom->tag[c_ind] << " c1: " << atom->tag[c1] << "\n";
				std::cout << "(dx, dy, dz): (" << dx << ", " << dy << ", " << dz << "\n"
					<< " Fdr = " << (fx*dx+fy*dy+fz*dz)/sqrt(dx*dx+dy*dy+dz*dz); 
				std::cout << " r[" << atom->tag[c_ind] << "]: " << x[c_ind][0] << ", " << x[c_ind][1] << ", " << x[c_ind][2] << ", "
					<<  " r[" << atom->tag[c1] << "]: " << x[c1][0] << ", " << x[c1][1] << ", " << x[c1][2] << "\n";
				*/
				fforce = sqrt(fx * fx + fy * fy + fz * fz);

				// dr = r_A0 dot r_10 is the distance along the bond of pol from c_ind
				double dr = dotsign[itype][jtype] * ( delx * x10 + dely * y10 + delz * z10);
				if( atom->tag[c_ind] > atom->tag[c1]){
					// Pulling towards c1 not c_ind, need to subtract from length of bond (|r_10|)
					dr = sqrt( x10 * x10 + y10 * y10 + z10 * z10) - dr;
				}

				double lsq = lx * lx + ly * ly + lz * lz;

				phi = fmag[itype][jtype] * dr + 0.5 * k_at[itype][jtype] * lsq - offset[itype][jtype];
				//std::cout << "Finished single compute\n";
				return factor_lj*phi;
			}else{
				fx=0.;
				fy=0.;
				fz=0.;
				fpair=0.;
				return 0.;
			}	    
		}
	}
	fx=0.;
	fy=0.;
	fz=0.;
	fpair=0.;
	return 0.;
}


