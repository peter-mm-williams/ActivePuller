import numpy as np
from lammps import lammps
import argparse
import random
from post_proc_funcs import *
import ctypes
from mpi4py import MPI

parser = argparse.ArgumentParser(description='Make a Polydisperse Packing')
parser.add_argument('-x', '--xyzfile', nargs='?',default='DNArelax.xyz', help='Output xyz file name', type=str)
parser.add_argument('-o', '--logfile', nargs='?',default='DNArelax.log', help='Output log file name', type=str)
parser.add_argument('-O', '--shortDat', nargs='?',default=0, help='Flag for if the simulation prints out short time data', type=int)
parser.add_argument('-m', '--molfile', nargs='?',default='DNA.mol', help='molecule file name', type=str)
parser.add_argument('-n', '--Ndt', nargs='?', default=5000, help='Number of dts per minimization', type=int)
parser.add_argument('-d', '--dt', nargs='?', default=0.0002, help='timestep', type=float)
parser.add_argument('-P', '--P0', nargs='?', default=1e-2, help='Desired Pressure', type=float)
parser.add_argument('-F', '--f0', nargs='?', default=0.5, help='Force of extrusion', type=float)
parser.add_argument('-f', '--phi', nargs='?', default=0.10, help='Initial Packing Fraction', type=float)
parser.add_argument('-T', '--T', nargs='?', default=1.0, help='Temperature', type=float)
parser.add_argument('-e', '--epsilon', nargs='?', default=1000.0, help='Well depth', type=float)
parser.add_argument('-s', '--sigCrowd', nargs='?', default=2.0, help='Diameter of Crowder particle in Pol complex', type=float)
parser.add_argument('-p', '--theta', nargs='?', default=132.0917709146271, help='Bond Angle', type=float)
parser.add_argument('-S', '--seed', nargs='?', default=1, help='Seed for pseudorandom number generator', type=int)
parser.add_argument('-t', '--Nsteps', nargs='?', default=2e7, help='Number of Steps for Simulation', type=int)
parser.add_argument('-l', '--lef', nargs='?', default=400, help='Scale for LEFs', type=int)
parser.add_argument('-b', '--bondangles', nargs='?', default=0, help='boolean for if system contains bond angles', type=int)
parser.add_argument('-k', '--k_att', nargs='?', default=5000.0, help='strength of Pol attraction', type=float)
parser.add_argument('-R', '--f_rand', nargs='?', default=0.0, help='Fraction of active force that is random', type=float)
parser.add_argument('-D', '--sig_p', nargs='?', default=0.6, help='range of Pol interaction', type=float)
parser.add_argument('-N', '--N', nargs='?', default=1, help='Number of Pols', type=int)
parser.add_argument('-X', '--si_lmp', nargs='?', default=0, help='Silence lammps output', type=int)
parser.add_argument('-B', '--BCs', nargs='?', default=0, help='Boundary Conditions: 0 for periodic cube, 1 for repulsive spherical', type=int)
parser.add_argument('-I', '--Int', nargs='?', default=1, help='Integration Schema: 0 for nve, 1 for overdamped langevin, 2 for langevin', type=int)

args = parser.parse_args()

def Pol_Sim(xyzfile, logfile, Ndt, dt, P_thresh, fmag, phi_target, seedval, Nsteps, T, epsilon, molfiles, theta, lef, BA, sigma3, k_att, sig_p, Npol, BCs, Integ, f_rand, shortDat, si_lmp):
	# Set up mpi4py

	import os
	from time import time
	t0=time()
	comm = MPI.COMM_WORLD
	me = comm.Get_rank()
	nprocs = comm.Get_size()
	phi0=0.01
	if me==0:
		print('Got mpi4py set up')

	# Delete files of the same names as our outputs should they exst:
	if me==0:
		try:
			f=open(xyzfile,'w')
			f.close()
		except:
			print('xyzfile does not exist yet')
		try:
			f=open(logfile,'w')
			f.close()
		except:
			print('logfile does not exist yet')
		try:
			f=open(xyzfile[:-3]+'dat','w')
			f.close()
		except:
			print('datfile does not exist yet')
		try:
			f=open(xyzfile[:-4]+'_short.dat','w')
			f.close()
		except:
			print('datfile does not exist yet')

	# Define Constants
	ntype = 6
	dt0=0.1
	dphi=0.01
	U_0=1e-10
	f0=5e-13
	U0=U_0
	U_rat=1e-10
	P_rat=1e-10
	ndim=3
	gamma = 0.01
	T0=1.0

	# Set up variable for lammps if si_lmp lammps outputs will be silenced
	if si_lmp:
		slt_screen=["-screen","none"]
		lmp = lammps(cmdargs=slt_screen)
	else:
		lmp = lammps()
	dt_step=10000
	logfile = xyzfile[:-4]+'.log'

	## Initialize Simulation, define units and dimensions
	lmp.command('units lj')
	lmp.command('dimension 3')
	lmp.command('newton on')
	lmp.command('boundary p p p')
	lmp.command('atom_style full')

	## Specify Interactions with lammps features
	lmp.command('variable rc equal 1.0')
	if epsilon==0:
		lmp.command('pair_style none')
	else:
		lmp.command('pair_style hybrid lj/rep/cut 1.0')
		lmp.command('pair_modify mix arithmetic')
	lmp.command('bond_style harmonic')
	if BA:
		lmp.command('angle_style cosine/squared')
		lmp.command('special_bonds angle yes')
	L=100
	Lx = L
	Ly = L
	Lz = L

	
	## Define Simulation Box
	lmp.command('region box block '+str(-Lx/2)+' '+str(Lx/2)+' '+str(-Ly/2)+' '+str(Ly/2)+' '+str(-Lz/2)+' '+str(Lz/2)+' units box')
	lmp.command('create_box '+str(ntype)+' box bond/types 2 angle/types 1 extra/bond/per/atom 2 extra/angle/per/atom 3 extra/special/per/atom 10')

	## Create Atoms

	# Parse file of moleculefile names
	molnames=[]
	print('got to molfiles')
	f=open(molfiles)
	for line in f:
		print(line)
		molnames.append(line.strip())

	# Load and create DNA molecule/ polymer
	for i in np.arange(0,len(molnames)):
		print(('molecule polymer'+str(int(i))+' '+molnames[i]))
		lmp.command('molecule polymer'+str(int(i))+' '+molnames[i])
	for i in np.arange(0,len(molnames)):
		lmp.command('create_atoms '+str(i)+' single 0.0 0.0 0.0 mol polymer'+str(int(i))+' '+str(seedval)+' rotate 0.0 1 0 0 ')
		print('create_atoms '+str(i)+' single 0.0 0.0 0.0 mol polymer'+str(int(i))+' '+str(seedval)+' rotate 0.0 1 0 0 \n')
	i=len(molnames)-1
	polID = 'polymer'+str(i)


	
	Natoms = lmp.get_natoms()
	nper_poly = Natoms/2
	Nmon = Natoms

	# Create vectors for use in post processing
	Nper_type=np.array([Natoms])
	masses=np.array([1.0])
	diameters=np.array([1.0])
	if me==0:
		print('\n\n\n THERE ARE '+str(Natoms)+' ATOMS IN THIS SYSTEM \n\n\n')


	## Define Masses and Pair Interactions
	sigma=1.0
	sig_arr = np.array([sigma, sig_p, sigma3, sigma, sigma, sigma, sig_p, sig_p, sig_p, sig_p, sig_p])
	masses = np.array([1,(sigma+sigma3)**3/2,(sigma+sigma3)**3/2,sigma**3, sigma**3, sigma**3])
	for i in np.arange(1,ntype+1):
		lmp.command('mass '+str(i)+' '+str(masses[i-1]))

	ind00 = 1
	ind01 = Nmon
	IndRange = np.array([ind00, ind01])
	
	for ind in np.arange(1,ind00+1):
		lmp.command('set atom '+str(ind)+' type 4')
	for ind in np.arange(ind01,Natoms+1):
		lmp.command('set atom '+str(ind)+' type 4')
	
	if epsilon!=0:
		for i in np.arange(1,ntype+1):
			for j in np.arange(i,ntype+1):
				lmp.command('pair_coeff '+str(i)+' '+str(j)+' lj/rep/cut '+str(epsilon)+' '+
					str((sig_arr[i-1]+sig_arr[j-1])/2.0)+' '+str((sig_arr[i-1]+sig_arr[j-1])/2.0))

	# Define Bond Interactions
	lmp.command('bond_coeff 1 '+str(3.0*T*10/0.08**2)+' 0.8 ')
	lmp.command('bond_coeff 2 '+str(3.0*T*10/((sigma3+1.0)/20.0)**2)+' '+str((sigma3+1.0)/2.0))
	if BA:
		lmp.command('angle_coeff 1 100.0 '+str(theta))
		lmp.command('special_bonds lj 0.0 0.0 1.0 angle yes dihedral no')
	else:
		lmp.command('delete_bonds all angle 1')
		lmp.command('special_bonds lj 0.0 1.0 1.0 angle no dihedral no')

	# Initialize Velocities and Neighborlist
	lmp.command('velocity all create '+str(T)+' '+str(seedval)+' mom no dist gaussian')
	#lmp.command('neighbor 0.8 bin')
	if epsilon==0:
		lmp.command('atom_modify sort 1000 0.5')


	# Relax System
	if BA:
		lmp.command('thermo_style custom step etotal epair ebond eangle pe ke fnorm fmax press')
	else:
		lmp.command('thermo_style custom step etotal epair ebond pe ke fnorm fmax press')

	lmp.command('thermo '+str(np.round(Ndt)))
	lmp.command('timestep '+str(dt))
	lmp.command('run 0')

	lmp.command('min_style fire')
	lmp.command('timestep '+str(dt0))
	lmp.command('minimize '+str(U0)+' '+str(f0)+' '+str(Ndt)+' '+str(Ndt))
	
	
	# Shrink box to desired packing fraction
	overlap = 0.8
	Volume = Natoms*np.pi/6.- 2*(nper_poly-1)*np.pi/12.*(2+overlap)*(1-overlap)**2 # Volume of polymer
	width = (Volume/phi_target)**(1./3)
	L=set_wall(Volume, phi_target ,lmp, 3)
	if BCs==0:
		width = (Volume/phi_target)**(1./3)
		print(Volume)
		print(phi_target)
		print(width)
		print(Volume/width**3)
	elif BCs==1:
		width = (1/phi_target*Volume*6./np.pi)**(1./3)/np.sqrt(2)
	if me==0:
		if BCs==1:
			print('width: %f, phi: %f, Volume: %f, phi_check: %f' %(width,phi_target, Volume, Volume/(np.pi*width**3/6.)))
		else:
			print('width: %f, phi: %f, Volume: %f, phi_check: %f' %(width,phi_target, Volume, Volume/(width**3)))

	lmp.command('minimize '+str(U0)+' '+str(f0)+' '+str(Ndt)+' '+str(Ndt))
	
	if BCs==1:
		lmp.command('region nucleus sphere 0.0 0.0 0.0 '+str((width)/2.0)+' units box side in')
		lmp.command('fix WallRep2 all wall/region nucleus lj126 1.0 '+str(1.0)+' '+str(2.**(1.0/6.0)))
	
	# Minimize structure twice
	lmp.command('timestep '+str(dt0))
	lmp.command('minimize '+str(U0)+' '+str(f0)+' '+str(Ndt)+' '+str(Ndt))
	lmp.command('minimize '+str(U0)+' '+str(f0)+' '+str(Ndt)+' '+str(Ndt))


	# ---------------------------- Equilibrate Chain under langevin thermostat ------------------------
	# Set up simulation for equilibrated chain
	lmp.command('fix 1 all nve')
	lmp.command('fix ln all langevin '+str(T0)+' '+str(T0)+' '+str(gamma*masses[int(0)])+' '+str(seedval))
	
	lmp.command('minimize '+str(U0)+' '+str(f0)+' '+str(Ndt)+' '+str(Ndt))
	
	# Set timestep
	lmp.command('timestep '+str(dt))

	# Unfix equilibrium values
	lmp.command('unfix 1')
	lmp.command('unfix ln')
	
	# Set up thermo for experiment simulation
	lmp.command('thermo '+str(np.round(Nsteps/1000)))
	
	# Change boundary conditions to shrinkwrapped
	lmp.command('change_box all boundary s s s')

	# ---------------------------- Add active puller molecules to system and set interactions ---------------------

	# Set polymerase interactions	
	set_pol_Interactions(lmp, fmag, sig_arr, k_att, ntype, Npol, epsilon, IndRange, f_rand, seedval)

	# Set up compute to access unboxed positions
	lmp.command('compute 1 all property/atom xu yu zu vx vy vz id proc')
	cval = lmp.extract_compute("1",1,2)

	# Make Pol molecules
	polatomIDs=np.array([])
	print(' Create polymerases')
	for pol_ind in np.arange(0,Npol):
		polatomIDs = create_pol(lmp, comm, polID, IndRange, polatomIDs, sig_arr, seedval, ntype, fmag, k_att, epsilon, cval, L)
	print(' PolatomIDs: ')
	print(polatomIDs)

	# Define pol interactions again after polymerases have been created
	set_pol_Interactions(lmp, fmag, sig_arr, k_att, ntype, Npol, epsilon, IndRange, f_rand, seedval)
	print('Defined interactions')
	
	# Define Groups
	for i in np.arange(1,ntype+1):
		grp_id='gr'+str(i)
		lmp.command('group '+grp_id+' type '+str(i))

	# Set up fixes for overdamped dynamics
	if Integ==0:
		lmp.command('fix 1 all nve')
	elif Integ==1:
		for i in np.arange(1,ntype+1):
			if i!=2 and i!=3:
				lmp.command('fix bd'+str(i)+' gr'+str(i)+' bd/baoab '+str(T)+' '+str(gamma*masses[int(i-1)])+' '+str(seedval))
			else:
				lmp.command('fix bd'+str(i)+' gr'+str(i)+' bd/baoab '+str(1e-4)+' '+str(gamma*masses[int(i-1)])+' '+str(seedval))
				print('fix bd'+str(i)+' gr'+str(i)+' bd/baoab '+str(1e-4)+' '+str(gamma*masses[int(i-1)])+' '+str(seedval)+'\n')
	else:
		lmp.command('fix 1 all nve')
		for i in np.arange(1, ntype+1):
			lmp.command('fix ln'+str(i)+' gr'+str(i)+' langevin '+str(T)+' '+str(gamma*masses[int(i-1)])+' '+str(seedval))

	Natoms = lmp.get_natoms()
		
	# Define frequency of printouts for log-spaced statistics
	l_0s = np.array([0.1,1.0,2.0])
	ks=2*np.pi/l_0s
	Nfrac=500
	Nper = Nsteps/Nfrac
	Nper_log = np.unique(np.round(np.logspace(0,np.log10(Nper),100)))
	Ndiffs = np.append(Nper_log[0], Nper_log[1:]-Nper_log[0:-1])
	
	totprint2=0
	lmp.command('compute ids all property/local patom1 patom2 ptype1 ptype2')
	lmp.command('compute fs all pair/local fx fy fz')
	lmp.command('timestep '+str(dt))
	
	t_ind=0
	t_check=dt_step
	v_pol=[]
	t1=time()

	# -------------------------------- Equilibrate System -----------------------------------------
	print('Equilibration:')

	lmp.command('dump 2 all custom '+str(int(np.round(1e7/100)))+' '+xyzfile[:-4]+'_equil.dat id type xu yu zu vz vy vz')
	lmp.command('run 10000000')
	lmp.command('undump 2')
	print('Done equilibrating')

	t2=time()
	
	if me==0:
		print('Run Experiment')
	# Set up dumps
	lmp.command('dump 1 all local '+str(int(np.round(Nsteps/Nfrac)))+' '+xyzfile[:-4]+'.pair index c_ids[1] c_ids[2] c_ids[3] c_ids[4] c_fs[1] c_fs[2] c_fs[3]')
	lmp.command('dump 2 all custom '+str(int(np.round(Nsteps/Nfrac)))+' '+xyzfile[:-4]+'.dat id type xu yu zu vz vy vz')
	lmp.command('dump_modify 2 sort id')
	
	if shortDat:
		lmp.command('dump 3 all local '+str(int(np.round(10*Nsteps)))+' '+xyzfile[:-4]+'_short.pair index c_ids[1] c_ids[2] c_ids[3] c_ids[4] c_fs[1] c_fs[2] c_fs[3]')
		lmp.command('dump_modify 3 first yes')
		lmp.command('write_dump all custom '+xyzfile[:-4]+'_short.dat id type xu yu zu vx vy vz modify sort id append yes')
		lmp.command('run 0')
	
	j=0
	
	for nt in np.arange(0,Nfrac):
		print('Done with %d out of %d steps' %(nt, Nfrac))
		if shortDat==1:
			for ind in np.arange(0,len(Ndiffs)):
				lmp.command('run '+str(int(np.round(Ndiffs[ind]))))
				lmp.command('write_dump all custom '+xyzfile[:-4]+'_short.dat id type xu yu zu vx vy vz modify sort id append yes')
				lmp.command('undump 3')
				lmp.command('dump 3 all local '+str(int(np.round(10*Nsteps)))+' '+xyzfile[:-4]+'_short.pair index c_ids[1] c_ids[2] c_ids[3] c_ids[4] c_fs[1] c_fs[2] c_fs[3]')
				lmp.command('dump_modify 3 first yes append yes')
		else:
			lmp.command('run '+str(Nper))
	if shortDat==1:
		lmp.command('undump 3')
		lmp.command('dump 3 all local '+str(int(Nsteps))+' '+xyzfile[:-4]+'_short.pair index c_ids[1] c_ids[2] c_ids[3] c_ids[4] c_fs[1] c_fs[2] c_fs[3]')
		lmp.command('dump_modify 3 first yes append yes')
		lmp.command('run 0')
	
	t3=time()
	if me==0:
		print('Finished Simulation, now processing. time of initialization: %f, time of equilibration: %f, time of simulation: %f, total time: %f' %((t1-t0),(t2-t1),(t3-t2),(t3-t0)))
	
	print("Proc %d out of %d procs has" % (me,nprocs), lmp)
	#MPI.Finalize()

def check_atBE(lmp, comm, IndRange, polatomIDs, cut, sig_arr, ntype, fmag, k_att, rs, L, epsilon, cval):
	print('In Check at BE')
	ndim=3
	for i in polatomIDs:
		print(int(i-1))
		dRs = rs[int(i-1),:] - rs[:100,:]
		#dRs = dRs - L*np.round(dRs/L)
		dR = np.sum((dRs)**2,axis=1)
		indmon = np.argmin(dR)
		#print('indmon: %d' %indmon)
		#print('Polymerase has coordinates: (%f, %f, %f) '%(rs[int(i-1),0],rs[int(i-1),1],rs[int(i-1),2]))
		#print(' Closest to %d which has coordinates: (%f, %f, %f)' %(indmon, rs[indmon,0],rs[indmon,1],rs[indmon,2]))
		#print(' Closest to %d which has coordinates: (%f, %f, %f)' %(indmon+1, rs[int(indmon+1),0],rs[int(indmon+1),1],rs[int(indmon+1),2]))
		#print(' Closest to %d which has coordinates: (%f, %f, %f)' %(indmon-1, rs[int(indmon-1),0],rs[int(indmon-1),1],rs[int(indmon-1),2]))
		#print(dR[indmon])
		if dR[indmon] > cut:
			print(' Polymerase with id: %d fell off' %(i))
			move_pol(lmp, IndRange, i, sig_arr, ntype, fmag, k_att, comm, len(polatomIDs), epsilon, lmp.get_natoms(), cval)
			print('Get new rs')
			lmp.command('run 1')
			rs, vs, N_cnt, which_proc = get_rsvs(cval, lmp.get_natoms(), ndim, comm, lmp)
			r_new=np.copy(rs[int(i-1),:])
			print('Check: new r_pol: '+str(r_new[0])+', '+str(r_new[1])+', '+str(r_new[2])+' ')
		elif indmon==np.min(IndRange) or indmon== np.max(IndRange):
			print(' Polymerase with id: %d is not on transcibing region of chain, local ind: %d' %(i,indmon))
			move_pol(lmp, IndRange, i, sig_arr, ntype, fmag, k_att, comm, len(polatomIDs), epsilon, lmp.get_natoms(), cval)
			print('Get new rs')
			lmp.command('run 1')
			rs, vs, N_cnt, which_proc = get_rsvs(cval, lmp.get_natoms(), ndim, comm, lmp)
			r_new=np.copy(rs[int(i-1),:])
			print('Check: new r_pol: '+str(r_new[0])+', '+str(r_new[1])+', '+str(r_new[2])+' ')
	print('Finished check at BE')

def set_wall(V, phi ,lmp, ndim):
	L=(V/phi)**(1.0/ndim)
	L_str=str(-L/2)+' '+str(L/2)
	change_box='change_box all x final '+L_str+' y final '+L_str+' z final '+L_str+' remap units box '
	lmp.command(change_box)
	return L

def pol_vs(polatomIDs, vs, rs, maxmon, L):
	v_pol = []
	for polatomID in polatomIDs:
		dr = get_drPolInd(maxmon, polatomID, rs, L)
		if not np.isnan(np.sum(dr)):
			v_pol.append(np.dot(vs[int(polatomID-1),:],dr)/np.linalg.norm(dr))
	return v_pol

def get_drPolInd(maxmon, polatomID, rs, L):
	dRs = rs[int(polatomID-1),:] - rs[:maxmon,:]
	dRs = dRs - L*np.round(dRs/L)
	dR = np.sum((dRs)**2,axis=1)
	indmon = np.argmin(dR)
	dr = np.zeros(3)
	if indmon>maxmon-2 or indmon<1:
		dr = np.nan
	elif dR[indmon+1]>dR[indmon-1]:
		dr = rs[indmon,:] - rs[indmon-1,:]
	else:
		dr = rs[indmon+1,:] - rs[indmon,:]
	return dr

def create_pol(lmp, comm, polID, IndRange, polatomIDs, sig_arr, seedval, ntype, fmag, k_att, epsilon, cval,Ls):
	Natoms = lmp.get_natoms()
	nextID = Natoms
	print('In create_pol')
	P, U0 = get_PU(lmp, comm)
	r_loc = np.array([0., 0., 0.])
	print(lmp.get_natoms())
	
	while lmp.get_natoms() == Natoms:
		r_gen = np.random.rand(3)*Ls - Ls/2
		lmp.command('create_atoms 2 single '+str(r_gen[0])+' '+str(r_gen[1])+' '+str(r_gen[0])+' units box')
		print('create_atoms 2 single '+str(r_gen[0])+' '+str(r_gen[1])+' '+str(r_gen[0])+' units box')
		print(lmp.get_natoms())
	
	print(lmp.get_natoms())

	polatomID = nextID+1
	polatomIDs = np.append(polatomIDs, np.array([polatomID]))
	print('Move pol')
	move_pol(lmp, IndRange, polatomID, sig_arr, ntype, fmag, k_att, comm, int(len(polatomIDs)+1), epsilon, Natoms, cval, randplace=True, U0=(U0+fmag))
	print('Moved')

	print('Made Pol')
	return polatomIDs

def move_pol(lmp, IndRange, polatomID, sig_arr, ntype, fmag, k_att, comm, Npol, epsilon, Natoms, cval, randplace=False, U0=np.nan):
	dE = 5
	print(' In move_pol')
	sigma3 = sig_arr[2]
	#print('got sigma3')
	# Half length of Pol Molecule
	dL = (1.0+sigma3)/2.0
	#if randplace:
	#	set_pol_Attraction(lmp, sig_arr, ntype, epsilon)
	mute_pol_interaction(lmp)
	#set_pol_Interactions_mute(lmp, fmag, sig_arr, k_att, ntype, Npol, epsilon)
	#print( ' mute pol interactions')
	#print( ' Get U')
	if randplace:
		lmp.command('run 0')
	P, U0 = get_PU(lmp,comm)
	#print('After run 0: U0: %f' %(U0))
	U=1e5
	ndim=3
	#print('Got U and P')
	#lmp.command('run 0')
	rs, vs, N_cnt, which_proc = get_rsvs(cval, lmp.get_natoms(), ndim, comm, lmp)
	#xs=lmp.gather_atoms("xu",1,ndim)
	#print(' Gathered coordinates ')

	print('Natoms: %d' %(lmp.get_natoms()))
	ind_pol= int(polatomID-1)
	ind_crowd = int(polatomID)
	r_pol0=np.copy(rs[ind_pol,:])
	bad_coords=True
	cnt=-1
	while bad_coords:
		cnt+=1
		#print(' Try to place '+str(polatomID)+' U0: '+str(U0)+', U: '+str(U))
		# Generate locus for placement
		if randplace:
			ind_loc = int(np.round(IndRange[0]-1+np.random.rand()*(IndRange[1]-IndRange[0])))
		else:
			ind_loc = int(np.round(IndRange[1]-2-1+np.random.rand()*3))
		ind_loc = comm.bcast(ind_loc, root=0)
		ind_next = ind_loc-1
		# Collect Atoms in list
		ndim=3
		r0 = np.copy(rs[int(ind_loc-1),:])
		r1 = np.copy(rs[int(ind_next-1),:])
		
		r_loc = (r0+r1)/2
		u_vec = r1 - r0
		u_vec /= np.linalg.norm(u_vec)
		u = np.random.rand(3)
		u /= np.linalg.norm(u)
		u_perp = np.cross(u_vec, u)
		u_perp /= np.linalg.norm(u_perp)
		u_perp = comm.bcast(u_perp, root=0)
		r_pol = np.copy(r_loc)
		r_crowd = r_pol + dL * u_perp

		lmp.command('set atom '+str(int(polatomID))+' x '+str(r_pol[0])+' y '+str(r_pol[1])+' z '+str(r_pol[2]))
		print('set atom '+str(int(polatomID))+' x '+str(r_pol[0])+' y '+str(r_pol[1])+' z '+str(r_pol[2]))
		if ~randplace:
			lmp.command('run 0')
		#print('ran 0')
		P, U = get_PU(lmp, comm)
		bad_coords = U>U0 + dE
		if cnt>100:
			bad_coords = False
			r_pol = np.copy(rs[int(IndRange[0]),:])
			lmp.command('set atom '+str(polatomID)+' x '+str(r_pol[0])+' y '+str(r_pol[1])+' z '+str(r_pol[2]))
			print('Give up set to where it originally was \n set atom '+str(polatomID)+' x '+str(r_pol[0])+' y '+str(r_pol[1])+' z '+str(r_pol[2]))
			ind_loc = IndRange[0]
			#lmp.command('run 0')
	rs, vs, N_cnt, which_proc = get_rsvs(cval, lmp.get_natoms(), ndim, comm, lmp)

	print('Placed: '+str(polatomID)+' U0: '+str(U0)+', U: '+str(U)+' in between '+str(ind_loc)+' and '+str(ind_next))
	print('r[%d]: %f, %f, %f' %(ind_loc,rs[int(ind_loc-1),0],rs[int(ind_loc-1),1],rs[int(ind_loc-1),2]))
	print('r[%d]: %f, %f, %f' %(ind_next,rs[int(ind_next-1),0],rs[int(ind_next-1),1],rs[int(ind_next-1),2]))
	print('r[%d]: %f, %f, %f' %(polatomID,rs[int(polatomID-1),0],rs[int(polatomID-1),1],rs[int(polatomID-1),2]))

	activate_pol_interaction(lmp)
	#set_pol_Interactions(lmp, fmag, sig_arr, k_att, ntype, Npol, epsilon)
	
	print(' r0: (%f, %f, %f), r1: (%f, %f, %f), r_loc: (%f, %f, %f)' %(r0[0],r0[1],r0[2], r1[0],r1[1],r1[2], r_pol[0],r_pol[1],r_pol[2]))
	P, U = get_PU(lmp, comm)
	print(' With Polymerase Interactions: U0: '+str(U0)+', U: '+str(U))

def set_pol_Interactions(lmp, fmag, sig_arr, k_att, ntype, Npol, epsilon, IndRange, f_rand, seedval):
	print('In Pol Interactions')
	lmp.command('pair_style hybrid lj/rep/cut '+str(1.0)+' POL 1.0')
	lmp.command('pair_coeff 1 1 lj/rep/cut '+str(epsilon)+' '+str(sig_arr[0])+' '+str(sig_arr[0]))
	print('pair_coeff 1 2 POL '+str(fmag)+' 1 '+str(k_att)+' '+str(Npol)+' '+str(int(IndRange[0]+1))+' '+
		str(IndRange[1])+' '+str(sig_arr[1])+' '+str(f_rand)+' '+str(3.*(sig_arr[0]+sig_arr[1])/2.0))
	lmp.command('pair_coeff 1 2 POL '+str(fmag)+' 1 '+str(k_att)+' '+str(Npol)+' '+str(int(IndRange[0]+1))+' '+
		str(IndRange[1])+' '+str(sig_arr[1])+' '+str(f_rand)+' '+str(seedval)+' '+str(3.*(sig_arr[0]+sig_arr[1])/2.0))
	lmp.command('pair_coeff 1 3 lj/rep/cut '+str(epsilon)+' '+str((sig_arr[0]+sig_arr[2])/2.0)+' '+str((sig_arr[0]+sig_arr[2])/2.0))
	lmp.command('pair_coeff 1 4 lj/rep/cut '+str(epsilon)+' '+str((sig_arr[0]+sig_arr[3])/2.0)+' '+str((sig_arr[0]+sig_arr[3])/2.0))
	lmp.command('pair_coeff 1 5 lj/rep/cut '+str(epsilon)+' '+str((sig_arr[0]+sig_arr[4])/2.0)+' '+str((sig_arr[0]+sig_arr[4])/2.0))
	lmp.command('pair_coeff 1 6 lj/rep/cut '+str(epsilon)+' '+str((sig_arr[0]+sig_arr[5])/2.0)+' '+str((sig_arr[0]+sig_arr[5])/2.0))

	print('Made all pair_coeffs with type 1')
	for i in np.arange(2,ntype+1):
		for j in np.arange(i,ntype+1):
			if i==3 and j==3:
				print('pair_coeff '+str(i)+' '+str(j)+' lj/rep/cut 0.001 '+
				str((sig_arr[i-1]+sig_arr[j-1])/2.0)+' '+str((sig_arr[i-1]+sig_arr[j-1])/2.0))
				lmp.command('pair_coeff '+str(i)+' '+str(j)+' lj/rep/cut 0.001 '+
					str((sig_arr[i-1]+sig_arr[j-1])/2.0)+' '+str((sig_arr[i-1]+sig_arr[j-1])/2.0))
			else:
				print('pair_coeff '+str(i)+' '+str(j)+' lj/rep/cut 1.0 '+
					str((sig_arr[i-1]+sig_arr[j-1])/2.0)+' '+str((sig_arr[i-1]+sig_arr[j-1])/2.0))
				lmp.command('pair_coeff '+str(i)+' '+str(j)+' lj/rep/cut 1.0 '+
					str((sig_arr[i-1]+sig_arr[j-1])/2.0)+' '+str((sig_arr[i-1]+sig_arr[j-1])/2.0))
	print('Set Pol Interactions')

def mute_pol_interaction(lmp):
	lmp.command('neigh_modify exclude type 1 2')

def activate_pol_interaction(lmp):
	lmp.command('neigh_modify exclude none')

def delete_molecule(lmp, atomID):
	'''
	Delete molecule of which atomID is the ID of a member atom
	'''
	lmp.command('group grID id '+str(atomID))
	lmp.command('delete_atoms group grID mol yes')
	lmp.command('group grID delete')

def get_forces(lmp, Natoms, comm, ndim):
	rank=comm.Get_rank()
	nprocs=comm.Get_size()

	fs=lmp.gather_atoms("f",1,ndim)
	f_max=0
	Nper=np.round(Natoms*ndim/nprocs)
	low_ind=int(rank*Nper)
	high_ind=int((rank+1)*Nper)
	if rank==nprocs-1:
		high_ind=Natoms*ndim
	for i in np.arange(low_ind,high_ind):
		if fs[i]>f_max:
			f_max=fs[i]
	f_maxes=comm.allgather(f_max)
	f_max=max(f_maxes)

	return f_max


def get_PU(lmp, comm):
	nprocs=comm.Get_size()
	U=lmp.extract_compute('thermo_pe',0,0)
	Us=comm.allgather(U)
	U=0
	for i in np.arange(0,nprocs):
		U+=Us[i]
	P=lmp.extract_compute('thermo_press',0,0)
	Ps=comm.allgather(P)
	P=0
	for i in np.arange(0,nprocs):
		P+=Ps[i]
	return P,U

def get_ribDist(Nrib, seedval, polysome_distfile):
    Dist=np.round(np.loadtxt(polysome_distfile))
    CumDist = np.cumsum(Dist/np.sum(Dist))

    np.random.seed(seedval)
    Nper = np.arange(1,8)
    Nmade_rib=0
    Nlist = np.zeros(len(Nper))
    while Nmade_rib < 0.63*Nrib:
        N_i = int(Nper[np.where(np.random.rand()<CumDist)[0][0]])
        if N_i>1:
            Nlist[N_i]+=1
            Nmade_rib+=N_i
    Nfree = Nrib - Nmade_rib
    return Nfree, Nlist

def get_rsvs(cval, Natoms, ndim, comm, lmp):
	#print('In get_rsvs')
	cval = lmp.extract_compute("1",1,2)
	rs = np.zeros((Natoms,ndim))
	vs = np.zeros((Natoms,ndim))
	me = comm.Get_rank()
	N_cnt = np.zeros((Natoms))
	which_proc = np.ones((Natoms))*1000
	for i in np.arange(0,lmp.extract_global("nlocal",0)):
		#print('Get ind')
		ind = int(cval[i][6]-1)
		#print('i: %d, ind: %d' %(i,ind))
		proc = int(cval[i][7])
		#print('Got proc data')
		if ind<0:
			print('Number of local atoms: %d' %lmp.extract_global("nlocal",0))
			print('Number of total atoms: %d' %lmp.get_natoms())
			print('ind is less than 0: %d, Exiting' %ind)
			exit(1)
		#if ind<0 or np.linalg.norm(rs[ind,:])>0:
		#	break
		#	print('me: %d, ind: %d, proc: %d' %(me,ind,proc))
		#print(proc)
		#if ind>=0 and me==proc:
		#	if np.round(cval[i][0], decimals=6) != 0.:
		#print(ind)
		N_cnt[ind]+=1
		which_proc[ind] = proc
		#print('filling rs: %d, with ind: %d' %(i,ind))
		rs[ind,:] = np.array([cval[i][0], cval[i][1], cval[i][2]])
		vs[ind,:] = np.array([cval[i][3], cval[i][4], cval[i][5]])
	return rs, vs, N_cnt, which_proc

def calc_theta_phi(r):
    '''
    Function that calculates the spherical polar angles (phi, psi) of a 
    vector in real space

    Inputs:
        r       position vector (vector of length 3 of floats)

    Outputs:
        theta   angle in radians defining declination from +z axis (float)
        phi     angle in radians from +x axis in xy plane (float)
    '''
    r_mag=np.linalg.norm(r)
    costheta=r[2]/r_mag
    theta=np.arccos(costheta)
    sintheta=np.sin(np.arccos(costheta))
    cosphi=r[0]/(r_mag*sintheta)
    sinphi=r[1]/(r_mag*sintheta)
    # Determine quadrant of angle based on signs of cos and sin and 
    #    place angle accordingly and convert to degrees
    if cosphi>0:
        phi=np.arcsin(sinphi)
    elif sinphi<0:
        phi=-np.arccos(cosphi)
    else:
        phi=np.arccos(cosphi)
    return theta, phi


Pol_Sim(args.xyzfile, args.logfile, args.Ndt, args.dt, args.P0, args.f0, args.phi, args.seed, args.Nsteps, args.T, args.epsilon, args.molfile, args.theta,
 args.lef, args.bondangles, args.sigCrowd, args.k_att, args.sig_p, args.N, args.BCs, args.Int, args.f_rand, args.shortDat, args.si_lmp)


