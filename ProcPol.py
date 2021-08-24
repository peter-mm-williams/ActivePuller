import numpy as np
from lammps import lammps
import argparse
import random
from post_proc_funcs import *
import ctypes

parser = argparse.ArgumentParser(description='Make a Polydisperse Packing')
parser.add_argument('-x', '--xyzfile', nargs='?',
                    default='DNArelax.xyz', help='Output xyz file name', type=str)
parser.add_argument('-o', '--logfile', nargs='?',
                    default='DNArelax.log', help='Output log file name', type=str)
parser.add_argument('-m', '--molfile', nargs='?',
                    default='DNA.mol', help='molecule file name', type=str)
parser.add_argument('-n', '--Ndt', nargs='?', default=5000,
                    help='Number of dts per minimization', type=int)
parser.add_argument('-d', '--dt', nargs='?', default=0.001,
                    help='timestep', type=float)
parser.add_argument('-P', '--P0', nargs='?', default=1e-2,
                    help='Desired Pressure', type=float)
parser.add_argument('-F', '--f0', nargs='?', default=0.5,
                    help='Force of extrusion', type=float)
parser.add_argument('-f', '--phi', nargs='?', default=0.10,
                    help='Initial Packing Fraction', type=float)
parser.add_argument('-T', '--T', nargs='?', default=1e-2,
                    help='Temperature', type=float)
parser.add_argument('-e', '--epsilon', nargs='?',
                    default=1.0, help='Well depth', type=float)
parser.add_argument('-s', '--sigCrowd', nargs='?', default=2.0,
                    help='Diameter of Crowder particle in Pol complex', type=float)
parser.add_argument('-p', '--theta', nargs='?',
                    default=141.871487167, help='Bond Angle', type=float)
parser.add_argument('-S', '--seed', nargs='?', default=1,
                    help='Seed for pseudorandom number generator', type=int)
parser.add_argument('-t', '--Nsteps', nargs='?', default=5e7,
                    help='Number of Steps for Simulation', type=int)
parser.add_argument('-l', '--lef', nargs='?', default=400,
                    help='Scale for LEFs', type=int)
parser.add_argument('-b', '--bondangles', nargs='?', default=0,
                    help='boolean for if system contains bond angles', type=int)
parser.add_argument('-k', '--k_att', nargs='?', default=50.0,
                    help='strength of Pol attraction', type=float)
parser.add_argument('-D', '--sig_p', nargs='?', default=0.6,
                    help='range of Pol interaction', type=float)
parser.add_argument('-N', '--N', nargs='?', default=1,
                    help='Number of Pols', type=int)
parser.add_argument('-M', '--Nmon', nargs='?', default=1,
                    help='Number of Monomers', type=int)
parser.add_argument('-B', '--BCs', nargs='?', default=0,
                    help='Boundary Conditions: 0 for periodic cube, 1 for repulsive spherical', type=int)

args = parser.parse_args()


def ProcPol(xyzfile, logfile, Ndt, dt, P0, fmag, phi_target, seedval, Nsteps, T, epsilon, molfiles, theta, lef, BA, sigma3, k_att, sig_p, Npol, BCs, Nmon):

    # Process .dat files

    # Load in Short time coordinate data from "_short.dat"
    Nfrac = 500
    filename = xyzfile[:-4]
    # breakpoint()
    # Get log spaced lag times
    Nper = Nsteps/Nfrac
    Nper_log = np.unique(np.round(np.logspace(0, np.log10(Nper), 100)))
    Ndiffs = np.append(Nper_log[0], Nper_log[1:]-Nper_log[0:-1])
    t, rs, vs, types, Ls, Natoms = process_lmp(xyzfile[:-4]+'_short.dat', dt)
    print('Processed short file')
    # breakpoint()

    # Calculate the Number of chain atoms in simulation
    Nchain = Natoms - Npol
    Nmon = Nchain
    if Npol > 0:
        x_pol, i1s1, i0s1, drijs1, drs1 = get_polpos(rs, Nchain, Nchain, Ls[0])
        tF, Fs = process_pair(filename+'_short.pair', 1e-3, Natoms, x_pol)
    print('Processed short pair file ')
    # Initialize variables for arrays of times for forces and coordinates
    ts = np.zeros((Nfrac, len(Ndiffs)))
    tsF = np.zeros((Nfrac, len(Ndiffs)))
    # breakpoint()

    # Fill in rs_master for positions of coordinates (Natoms x ndim x Number of data points per lag x Number of lags)
    rs_master = np.zeros((Natoms, 3, Nfrac, len(Ndiffs)))*np.nan
    i = 0
    for nt in np.arange(0, Nfrac):
        for ind in np.arange(0, len(Ndiffs)):
            ts[nt, ind] = t[i]
            rs_master[:, :, nt, ind] = rs[:, :, i]
            i += 1
            if i > np.shape(rs)[2]:
                break
        if i > np.shape(rs)[2]:
            break
    # save rs_master to npy file
    np.save(xyzfile[:-4] + '.Rshort.npy', rs_master)
    # breakpoint()

    # Load in forces at short times and calculate force correlation functions
    if Npol > 0:
        fs_master = np.zeros((Natoms, 3, Nfrac, len(Ndiffs)))
        i = 0
        for nt in np.arange(0, Nfrac):
            for ind in np.arange(0, len(Ndiffs)):
                tsF[nt, ind] = tF[i]
                fs_master[:, :, nt, ind] = Fs[:, :, i]
                i += 1
                if i > np.shape(Fs)[2]:
                    break
            if i > np.shape(Fs)[2]:
                break
        np.save(xyzfile[:-4]+'.Fshort.npy', fs_master)
        dis = np.arange(1, 5)
        FAC = np.zeros(len(Ndiffs))
        FPC = np.zeros((len(dis), 2*len(Ndiffs)))
        fmean = np.nanmean(np.nanmean(
            fs_master[:Nchain, :, :, 0], axis=2), axis=1)
        fmean2 = np.nansum(fmean*fmean)
        FACnorm = np.nanmean(
            np.sum(fs_master[:Nchain, :, :, 0]*fs_master[:Nchain, :, :, 0], axis=1))-fmean2
        for i in np.arange(0, len(Ndiffs)):
            FAC[i] = (np.nanmean(np.sum(fs_master[:Nchain, :, :, i]
                                        * fs_master[:Nchain, :, :, 0], axis=1))-fmean2)/FACnorm
        for i in np.arange(-len(Ndiffs), len(Ndiffs)):
            for j in np.arange(0, len(dis)):
                if i < 0:
                    FPC[j, i] = np.nanmean(np.sum(
                        fs_master[dis[j]:Nchain, :, :, i]*fs_master[:Nchain-dis[j], :, :, 0], axis=1))
                    f_norm = np.nanmean(
                        np.sqrt(np.sum(fs_master[dis[j]:Nchain, :, :, :]**2, axis=1)))
                    g_norm = np.nanmean(
                        np.sqrt(np.sum(fs_master[:Nchain-dis[j], :, :, :]**2, axis=1)))
                    FPC[j, i] /= (f_norm*g_norm)
                    # np.sqrt(np.mean(np.sum(fs_master[dis[j]:Nchain,:,:,:]**2,axis=1))*np.mean(np.sum(fs_master[:Nchain-1-dis[j],:,:,:]**2,axis=1)))
                else:
                    FPC[j, i] = np.nanmean(np.sum(
                        fs_master[dis[j]:Nchain, :, :, 0]*fs_master[:Nchain-dis[j], :, :, i], axis=1))
                    f_norm = np.nanmean(
                        np.sqrt(np.sum(fs_master[dis[j]:Nchain, :, :, :]**2, axis=1)))
                    g_norm = np.nanmean(
                        np.sqrt(np.sum(fs_master[:Nchain-dis[j], :, :, :]**2, axis=1)))
                    FPC[j, i] /= (f_norm*g_norm)
        np.save(xyzfile[:-4]+'.FPCshort.npy', FPC)
        np.save(xyzfile[:-4] + '.FACshort.npy', FAC)
        # breakpoint()
    # breakpoint()

    # Calculate correlation of end to end vector from displacement vector
    # initialize dR for array of displacements for short time data
    dR = np.copy(rs_master)*np.nan
    Ree = rs_master[Nmon, :, :, :]-rs_master[0, :, :, :]
    Ct = np.zeros((np.shape(rs_master)[2], np.shape(rs_master)[3]))*np.nan
    dts = np.zeros(np.shape(rs_master)[3])
    # Fill displacement array and end-to-end vector correlation function as a function of time lag
    for i in np.arange(0, np.shape(rs_master)[3]):
        dR[:, :, :, i] = rs_master[:, :, :, i]-rs_master[:, :, :, 0]
        dts[i] = ts[0, i]-ts[0, 0]
        Ct[:, i] = np.nanmean(Ree[:, :, i]*Ree[:, :, 0], axis=0)
    # Calculate mean square displacement
    MSD = np.nanmean(np.sum(dR**2, axis=1), axis=1)
    # Take mean and normalize end-to-end vector correlation function ("Ct_norm")
    Ct_norm = np.nanmean(Ct, axis=0)
    Ct_norm /= Ct_norm[0]
    np.save(xyzfile[:-4]+'.Ctshort.npy', Ct_norm)
    np.save(xyzfile[:-4]+'.MSDshort.npy', MSD)

    '''
	------------------------------------ Make calculations for long time data (.dat file) -----------------------------
	'''
    n_steps = np.unique(np.round(np.logspace(0, np.log10(Nsteps/500), 10)))
    Ndt = 400  # Goal number of logarithmically spaced lag times to calculate
    lmpfile = xyzfile[:-4]+'.dat'
    datfile = xyzfile[:-4]+'.txt'

    # Extract info from .dat file and save coordinates as rs.npy and end-to-end vector as R.npy
    t, rs, vs, types, Ls, Natoms = process_lmp(lmpfile, dt)
    np.save(xyzfile[:-4]+'.rs.npy', rs)
    np.save(xyzfile[:-4]+'.R.npy', np.array(rs[Nmon-1, :, :]-rs[0, :, :]))
    Ree = rs[Nmon-1, :, :]-rs[0, :, :]

    # Make .xyz file for visualization
    make_tclfile(datfile[:-4]+'.tcl', np.ones(len(np.unique(types))+1), Ls[0])
    write_xyzfile(datfile[:-4]+'.xyz', rs, types, Ls[0])

    # Generate logarithmically spaced indices for correlation functions
    inds = np.append(np.array([0]), np.unique(
        np.round(np.logspace(np.log10(1), np.log10(len(t)-2), Ndt)))).astype(int)
    print('inds:')
    print(inds)
    # Calculate correlation function of the end-to-end vector
    C_t = get_normalize_correlation_function(Ree.T, inds)
    np.save(xyzfile[:-4]+'.Ct.npy', C_t)
    np.save(xyzfile[:-4]+'.Cdt.npy', inds*(t[1]-t[0]))
    np.save(xyzfile[:-4]+'.t.npy', t)

    # Calculate radius of gyration of the chain and save to .Rg.npy
    N = np.shape(rs)[0]
    r_com = np.mean(rs[0:Nmon, :, :], axis=0)
    r_com = np.tile(r_com, (Nmon, 1, 1))
    Rg = np.sqrt(np.mean(np.sum((rs[:Nmon, :, :]-r_com)**2, axis=1), axis=0))
    np.save(xyzfile[:-4]+'.Rg.npy', Rg)

    # Calculate the radius of gyration for subchains of all lengths and save the mean and standard deviation
    Rgs_m = np.zeros((Nmon, 2))*np.nan
    print('Nmon: %d' % (Nmon))
    for i in np.arange(0, Nmon):
        print('i: %d, Nmon-i: %d' % (i, Nmon-i))
        Rgs_temp = np.array([])

        for j in np.arange(0, Nmon-i):
            #print('j: %d' %j)
            r_com = np.mean(rs[j:j+i, :, :], axis=0)
            Rg_temp = np.sqrt(
                np.mean(np.sum((rs[j:j+i, :, :]-r_com)**2, axis=1), axis=0))
            # print(Rg_temp)
            Rgs_temp = np.append(Rgs_temp, Rg_temp)
        print('Rg: %f' % np.nanmean(Rgs_temp))
        Rgs_m[i, 0] = np.nanmean(Rgs_temp)
        Rgs_m[i, 1] = np.nanstd(Rgs_temp)
    np.save(xyzfile[:-4]+'.Rgsm.npy', Rgs_m)

    # Calculate a histogram of the radius of gyration
    binvals = np.arange(0, Nmon, 0.25)
    Rghist, bin_edges = np.histogram(Rg[np.where(~np.isnan(Rg))], bins=binvals)
    bin_centers = 0.5 * (bin_edges[1:]+bin_edges[:-1])
    np.save(xyzfile[:-4]+'.Rghist.npy', Rghist)
    np.save(xyzfile[:-4]+'.Rgbin.npy', bin_centers)

    # Calculate the mean square displacement and the velocity autocorrelation function of the center of mass of the polymer
    inds = np.unique(np.round(np.logspace(
        np.log10(1), np.log10(len(t)-2), Ndt)))
    dts = inds*float(t[1]-t[0])
    MSDcom = np.zeros(len(inds))
    MFDcom = np.zeros(len(inds))
    VACcom = np.zeros(len(inds))
    vsCOM = np.nanmean(vs, axis=0)
    rsCOM = np.nanmean(rs, axis=0)
    for i in np.arange(0, len(inds)):
        ind = int(inds[i])
        drs = rsCOM[:Nmon, 0:-ind]-rsCOM[:Nmon, ind:]
        MSDcom[i] = np.mean(np.sum(drs**2, axis=0), axis=0)
        MFDcom[i] = np.mean(np.sum(drs**4, axis=0), axis=0)
        vdv = vsCOM[:Nmon, 0:-ind]*vsCOM[:Nmon, ind:]
        VACcom[i] = np.mean(np.sum(vdv, axis=0), axis=0)
    np.save(xyzfile[:-4]+'.vacCOM.npy', VACcom)
    np.save(xyzfile[:-4]+'.msdCOM.npy', MSDcom)
    np.save(xyzfile[:-4]+'.dtsCOM.npy', dts)

    # Calculate velocity autocorrelation function for each atom in system and save
    dts, VAC = make_VAC(t, vs, int(Nmon+Npol), Ndt)
    np.save(xyzfile[:-4]+'.vac_extra.npy', VAC)

    # Calculate the mean square displacement and the mean fourth moment of the displacement vector
    dts, MSD, MFD = make_MSDMFD(t, rs, np.shape(rs)[0], Ndt)
    np.save(xyzfile[:-4]+'.dt.npy', dts)
    np.save(xyzfile[:-4]+'.msd.npy', MSD)
    np.save(xyzfile[:-4]+'.mfd.npy', MFD)

    # Calculate the mean square displacement with the center of mass subtracted off
    dts, MSDsub, MFDsub = make_MSDMFD(t, rs[:Nmon, :, :]-r_com, Nmon, Ndt)
    np.save(xyzfile[:-4]+'.msdCOMsub.npy', MSDsub)
    np.save(xyzfile[:-4]+'.mfdCOMsub.npy', MSDsub)

    # Output costhetas.npy to calculate decay of chain (serves as measure for experimental persistence length)
    cos_thetas = get_costheta_ij(rs, Nmon)
    np.save(xyzfile[:-4]+'costhetas.npy', cos_thetas)

    # Make contact map from simulation sampling every 10 printout snapshots
    Nevery = 10
    for thresh in np.arange(1.0, 2.51, 0.5):
        HiC = make_HiC(rs, thresh, Nevery)
        np.save(xyzfile[:-4]+'.hic.'+str(thresh)+'.npy', HiC)

    # Get coarser image of system by averaging together 5 adjacent bins to singular bins
    N = Nmon
    Nevery = 1
    res = 5
    for thresh in np.arange(1.0, 3.0, 0.5):
        HiC = make_HiC(rs, thresh, Nevery)
        print('made HiC')
        HiCtrim = HiC[np.arange(0, N, res), :]
        HiCtrim = HiCtrim[:, np.arange(0, N, res)]
        for i in np.arange(0, np.shape(HiCtrim)[0]):
            for j in np.arange(0, np.shape(HiCtrim)[1]):
                HiCtrim[i, j] = np.nanmean(
                    HiC[i*res:(i+1)*res, j*res:(j+1)*res])
        np.save(xyzfile[:-4]+'_'+str(thresh)+'.hic.npy', HiC)
        np.save(xyzfile[:-4]+'_'+str(thresh)+'.hictrim.npy', HiCtrim)
    print('saved files')


def calc_lp(rs, b):
    '''
    Inputs:
            rs     position coordinates (N x ndim x Nt)
            b      scalar, relaxed length of bond of polymer
    Outputs:
            lp     persistence length as calulated by -b/<cos(theta)> where theta is the angle formed by 3 consecutive monomers in the polymer
    '''
    dr = rs[1:, :, :]-rs[:-1, :, :]
    costheta = np.sum(dr[:-1, :, :]*dr[1:, :, :], axis=1)/(
        np.linalg.norm(dr[:-1, :, :], axis=1)*np.linalg.norm(dr[1:, :, :], axis=1))
    lp = -b/np.mean(costheta, axis=0)
    return lp


def calc_Rg(rs):
    '''
    Function to calculate the radius of gyration

    Inputs:
            rs     position coordinates (N x ndim x Nt)
    Outputs:
            Rg     vector of the radius of gyration (length Nt)

    '''
    N = np.shape(rs)[0]
    r_com = np.tile(np.nanmean(rs, axis=0), (N, 1, 1))
    Rg = np.sqrt(1/N*np.sum(np.sum((rs-r_com)**2, axis=1), axis=0))
    return Rg


def get_HiCslice(r_s, thresh):
    '''
    Calculate a singulat time slice of a contact map
    Inputs:
            r_s       Array of position coordinates (N x ndim)
            thresh    Maximum distance to count as a contact
    Outputs:
            counts    N x N binary matrix of whether monomer i is in contact with monomer j for this configuration
    '''
    # print(np.shape(r_s))
    dx = np.subtract.outer(r_s[:, 0], r_s[:, 0])
    print('get dx')
    dy = np.subtract.outer(r_s[:, 1], r_s[:, 1])
    dz = np.subtract.outer(r_s[:, 2], r_s[:, 2])
    print('calc dists')
    dists = np.sqrt(dx**2+dy**2+dz**2)
    counts = dists < thresh
    return counts


def get_normalize_correlation_function(A, inds):
    '''
    Function to calculate the normalized autocorrelation function of variable A

    Inputs:
            A        array of a measure of the system (Nt x ndim)
            inds     vector of index lags to calculate correlation about
    Outputs:
            C        normalized autocorrelation vector (length number of lag times i.e. length(inds))
    '''
    C = np.zeros(len(inds))
    A2mean = np.mean(np.sum(A**2, axis=1))
    i = 0
    for ind in inds:
        if ind > 0:
            C[i] = np.mean(np.sum(A[ind:, :]*A[:-ind, :], axis=1))
        else:
            C[i] = np.mean(np.sum(A**2, axis=1))
    C /= A2mean
    return C


def process_pair(lmpfile, dt, Natoms, x_pol):
    '''
    Function to extract Forces from .pair file

    Inputs:
            lmpfile       path and filename for .pair file (string)
            dt            timestep of simulation (scalar float)
            Natoms        Number of atoms in simulation (scalar int)
            x_pol
    Outputs:
            t             vector of times (length Nt, float)
            Fs            array of forces
    '''
    x_lo = np.floor(x_pol).astype(int)
    x_hi = np.ceil(x_pol).astype(int)
    t = []
    Fs = []
    switch = 0
    Nval = 0
    cnt = 0
    tval = 0
    l_switch = 0
    ind = -1
    f = open(lmpfile)
    for x in f:
        if x.strip() == 'ITEM: TIMESTEP':
            ind += 1
            # t.append(float(f.readline().strip()))
            # print(ind)
            # print(t[-1])
            tval = 1
            pos_data = 0
            if switch > 0:  # and cnt%10==0:
                Fnorm = np.linalg.norm(fs, axis=1)
                # print(Fnorm)
                Fs = emptydstack(Fs, np.copy(fs))
                # print(fs)
                if switch == 1:
                    switch = 2
                fs = np.zeros((Natoms, 3))
        elif tval == 1:
            tval = 0
            # print(x.strip())
            t.append(float(x.strip()))
            # print(t[-1])
        elif x.strip() == 'ITEM: NUMBER OF ENTRIES' and switch == 0:
            #print('found number of entries')
            Nval = 1
        elif Nval == 1:
            Nval = 0
            N = int(x.strip())
        elif x[0:16] == 'ITEM: BOX BOUNDS' and switch == 0:
            l_switch = 1
        elif l_switch == 1:
            l_switch = 2
            l_line = x.strip().split(' ')
            # print(l_line)
            Lx = -float(l_line[0])+float(l_line[1])
            # print(float(l_line[0]))

            # print(Lx)
        elif l_switch == 2:
            l_switch = 3
            l_line = x.strip().split(' ')
            # print(l_line)
            Ly = -float(l_line[0])+float(l_line[1])
        elif l_switch == 3:
            l_switch = 0
            l_line = x.strip().split(' ')
            # print(l_line)
            Lz = -float(l_line[0])+float(l_line[1])
            Ls = np.array([Lx, Ly, Lz])
        elif x[0:13] == 'ITEM: ENTRIES':
            #print('found entries')
            cnt += 1
            pos_data = 1
            if switch == 0:
                switch = 1
                # print('initialize')
                # rs=np.zeros((N,3,1000))
                # vs=np.zeros((N,3,1000))
                #print('initialization complete')
                fs = np.zeros((Natoms, 3))
        elif pos_data == 1:
            line = x.strip().split(' ')
            atomid0 = int(line[1])-1
            atomid1 = int(line[2])-1
            atomtype0 = int(line[3])
            atomtype1 = int(line[4])

            fx = float(line[5])
            fy = float(line[6])
            fz = float(line[7])
            if (atomtype0 == 2 and atomtype1 == 1) or (atomtype1 == 2 and atomtype0 == 1):
                # print('before:')
                # print(fs)
                #print('x_lo: %d, x_hi:  %d' %(x_lo[ind],x_hi[ind]))
                #print('types: %d %d' %(atomtype0, atomtype1))
                #print('ids: %d %d' %(atomid0, atomid1))
                # print(np.array([-fx,-fy,-fz]))
                if atomtype0 == 2:
                    if atomid1 == x_lo[ind] or atomid1 == x_hi[ind]:
                        fs[atomid1, :] += np.array([-fx, -fy, -fz])
                        fs[atomid0, :] += np.array([fx, fy, fz])
                else:
                    if atomid0 == x_lo[ind] or atomid0 == x_hi[ind]:
                        fs[atomid0, :] += np.array([fx, fy, fz])
                        fs[atomid1, :] += np.array([-fx, -fy, -fz])
                # print('after:')
                # print(fs)
        # if ind>5:
        #    break
    Fnorm = np.linalg.norm(fs, axis=1)
    Fs = emptydstack(Fs, fs)
    t = np.array(t)*float(dt)
    print(t)
    return t, Fs


def get_polpos(rs, pol_ind, Nchain, L):
    '''
    Inputs:
            rs         array of coordinates (shape: N x ndim x Nt)
            pol_ind    index of puller (scalar int)
            Nchain     Number of monomers in polymer chain (scalar int)
            L          length of box edge (scalar float)
    Outputs:
            x_pol      vector for position of active puller along chain (length Nt, floats)
            i1s        higher index of monomer composing line segment nearest puller (length Nt, ints)
            i0s        lower index of monomer composing line segment nearest puller (length Nt, ints)
            drijs      displacement vectors pointing from higher index (i1s) to the active puller (pol_ind) (ndim x Nt, floats)
            drs        displacement vectors pointing from higher index (i1s) to lower index (i0s) (ndim x Nt, floats)
    '''
    Nt = np.shape(rs)[2]
    x_pol = []
    i1s = []
    i0s = []
    drijs = []
    drs = []
    for i in np.arange(0, Nt):
        dRvec = rs[:Nchain, :, i] - rs[pol_ind, :, i]
        #print('shape of dRvec')
        # print(np.shape(dRvec))
        #print('Nchain: %d, pol_ind %d' %(Nchain,pol_ind))
        dRvec -= L*np.round(dRvec/L)
        dR = np.linalg.norm(dRvec, axis=1)
        #print('shape of dR')
        # print(np.shape(dR))
        i_min = np.argmin(dR[1:Nchain-1])+1
        if int(i_min+1) >= Nchain:
            i1 = i_min
            i0 = i_min-1
        elif int(i_min) == 0:
            i1 = 1
            i0 = 0
        elif dR[int(i_min+1)] <= dR[int(i_min-1)]:
            i1 = i_min + 1
            i0 = i_min
            i1 = i0+1
        else:
            i1 = i_min
            i0 = i_min-1
        i1s.append(i1)
        i0s.append(i0)
        # Displacement vector pointing from higher index to lower index
        dr = rs[i0, :, i] - rs[i1, :, i]
        dr -= L*np.round(dr/L)
        # displacement vector pointing from higher index to pol
        drij = rs[pol_ind, :, i] - rs[i1, :, i]
        drij -= L*np.round(drij/L)
        drijs = emptydstack(drijs, drij)
        if i1 - np.dot(dr, drij)/np.linalg.norm(dr) < 0:
            print('i1 %d, i0: %d' % (i1, i0))
            print('Distance to high index: %f,  Distance to low index: %f' %
                  (dR[i1], dR[i0]))
            print('Distance to high index+1: %f' % (dR[i1+1]))
            u = dr/np.linalg.norm(dr)
            print('u: %f, %f, %f,  u.drij: %f' %
                  (u[0], u[1], u[2], np.dot(u, drij)))
            print('dR[min]')
            print(dR[i_min])
            print('x_pol: %f' % (i1 - np.dot(u, drij)))
            print(dR)
            break
        x_pol.append(i1 - np.dot(dr, drij)/np.linalg.norm(dr))
        drs = emptydstack(drs, dr)
    return np.array(x_pol), np.array(i1s), np.array(i0s), drijs, drs


def get_costheta_ij(rs, Nmon):
    '''
    Function that gets the cos of the angle of bond vector an index distance j away along the chain

    Inputs:
            rs            array of coordinates (shape: N x ndim x Nt)
            Nmon          Number of monomers
    Outputs:
            cos_thetas   vector of <cos(theta)> where theta = angle between bond vector i and vector i+j (length Nt of floats)
    '''

    Nt = np.shape(rs)[2]
    rij = rs[1:Nmon, :, :]-rs[:Nmon-1, :, :]
    rij_hat = np.copy(rij)
    cos_thetas = np.zeros(((Nmon-2), Nt))
    rij_norm = np.linalg.norm(rij, axis=1)
    for d in np.arange(0, 3):
        rij_hat[:, d, :] = rij[:, d, :]/rij_norm
    for ij in np.arange(0, Nmon-2):

        if ij > 0:
            print(np.shape(
                np.mean(np.sum(rij_hat[:-ij, :, :]*rij_hat[ij:, :, :], axis=1), axis=0)))
            print(np.shape(cos_thetas))
            print(np.shape(cos_thetas[ij, :]))
            cos_thetas[ij, :] = np.mean(
                np.sum(rij_hat[:-ij, :, :]*rij_hat[ij:, :, :], axis=1), axis=0)
        else:
            cos_thetas[ij, :] = np.mean(
                np.sum(rij_hat[:, :, :]*rij_hat[:, :, :], axis=1), axis=0)
    return cos_thetas


def make_HiC(rs, thresh, Nevery):
    '''
    Construct a contact map from rs of the coordinates in rs

    Input:
            rs        array of coordinates (shape: N x ndim x Nt)
            thresh    Maximum distance to count as a contact (scalar flt)
            Nevery    the number of time snapshots between sampling (scalar int)
    Outputs:
            counts    the raw counts of contacts between atoms in system(NxN of ints)
    '''
    Nsteps = np.shape(rs)[2]
    Natoms = np.shape(rs)[0]
    ind = 0
    counts = np.zeros((Natoms, Natoms))
    print('initialized counts')
    for i in np.arange(0, Nsteps, Nevery):
        ind += 1
        counts += get_HiCslice(rs[:, :, i], thresh)
        print(i)
    counts /= ind
    return counts


def get_FCC_FAC(fs, inds, dis):
    '''
    Extract cross correlation and autocorrelation of forces

    Inputs:
            fs     array of forces (N x ndim x Nt)
            inds   time index lags for autocorrelation (vector of ints)
            dis    atom index lags for cross correlation (vector of ints)
    Outputs:
            FCC    array of normalized cross correlation function for the force ()
            FAC    array of normalized auto correlation function for the force
    '''
    taus = inds.astype(int)
    dis = dis.astype(int)
    FCC = np.zeros((len(taus), len(dis)))
    j = -1
    for di in dis:
        j += 1
        fpc = np.zeros(len(taus))
        i = -1
        for tau in taus:
            i += 1
            fpc[i] = np.nanmean(
                np.mean(np.sum(fs[di:, :, tau:]*fs[:-di, :, :-tau], axis=1), axis=1))
        fpc /= np.sqrt(np.nanmean(np.sum(fs[di:, :, :] **
                                         2, axis=1)*np.sum(fs[:-di, :, :]**2, axis=1)))
        FCC[:, j] = fpc

    fac = np.zeros(len(taus))
    i = -1
    for tau in taus:
        i += 1
        fac[i] = np.nanmean(
            np.mean(np.sum(fs[:, :, tau:]*fs[:, :, :-tau], axis=1), axis=1))
    fac /= np.nanmean(np.sum(fs**2, axis=1))
    FAC = fac
    return FCC, FAC


def get_dr(rs, PolInd, Nmon, thresh, tind, Ls):
    '''
    Inputs:
            rs
            PolInd
            Nmon
            thresh
            tind
            Ls
    Outputs:
            dr
            mol0
    '''
    dRs = rs[PolInd, :, tind] - rs[:Nmon, :, tind]
    dRs -= Ls[0]*np.round(dRs/Ls[0])
    dR = np.linalg.norm(dRs, axis=1)
    mol0 = np.argmin(dR)
    print('mol0: %d, dR[mol0]: %f' % (mol0, dR[mol0]))
    # print(dR[mol0])
    if dR[mol0] < thresh and mol0 < Nmon-1:
        if dR[mol0+1] > dR[mol0-1]:
            # dr goes between mol0 and mol0-1
            dr = rs[mol0, :]-rs[mol0-1, :]
        else:
            dr = rs[mol0+1, :]-rs[mol0, :]
    else:
        dr = np.ones(3)*np.nan
    if dR[mol0] > thresh:
        mol0 = np.nan
    return dr, mol0


def get_dRt2(rs, mol0s, tind, taus, inds):
    dRt = np.zeros((len(inds), len(taus)))*np.nan
    i_t = -1
    for tau in taus:
        i_t += 1
        i = -1
        for ind in inds:
            i += 1
            # dRt[i,i_t]=np.linalg.norm(rs[int(mol0s[i_t]-ind),:,int(tind+tau)]-rs[int(mol0s[i_t]-ind),:,int(tind)])
            # print(np.linalg.norm(rs[int(mol0s[i_t]-ind),:,int(tind+tau)]-rs[int(mol0s[i_t]-ind),:,int(tind)]))
            try:
                dRt[i, i_t] = np.linalg.norm(
                    rs[int(mol0s[i_t]-ind), :, int(tind+tau)]-rs[int(mol0s[i_t]-ind), :, int(tind)])
            except:
                dRt[i, i_t] = np.nan
    return dRt


def get_dRtVec(rs, mol0s, tind, taus, inds):
    dRt0 = np.zeros((len(inds), len(taus), 3))*np.nan
    dRAC = np.zeros((len(taus)))
    dRCC = np.zeros((len(inds), len(taus)))
    dRt = np.zeros((len(inds), len(taus), 3))*np.nan
    i_t = -1
    for tau in taus:
        i_t += 1
        dR0 = rs[int(mol0s[i_t]-ind)]
        i = -1
        for ind in inds:
            i += 1
            # dRt[i,i_t]=np.linalg.norm(rs[int(mol0s[i_t]-ind),:,int(tind+tau)]-rs[int(mol0s[i_t]-ind),:,int(tind)])
            # print(np.linalg.norm(rs[int(mol0s[i_t]-ind),:,int(tind+tau)]-rs[int(mol0s[i_t]-ind),:,int(tind)]))
            try:
                dRt[i, i_t, :] = rs[int(
                    mol0s[i_t]-ind), :, int(tind+tau)]-rs[int(mol0s[i_t]-ind), :, int(tind)]
            except:
                dRt[i, i_t] = np.nan
    return dRt


def get_dRt(rs, mol0, tind, taus, inds):
    dRt = np.zeros((len(inds), len(taus)))*np.nan
    i_t = -1
    for tau in taus:
        i_t += 1
        i = -1
        for ind in inds:
            i += 1
            # print(np.linalg.norm(rs[int(mol0-ind),:,int(tind+tau)]-rs[int(mol0-ind),:,int(tind)]))
            try:
                dRt[i, i_t] = np.linalg.norm(
                    rs[int(mol0-ind), :, int(tind+tau)]-rs[int(mol0-ind), :, int(tind)])
            except:
                dRt[i, i_t] = np.nan
    return dRt


def get_dRts(rs, PolInd, Nmon, thresh, taus, inds, Ls):
    Nt = np.shape(rs)[2]
    max_t = Nt - np.max(taus)-1
    print(max_t)
    print(len(inds))
    print(len(taus))
    dRts = np.zeros((max_t, len(inds), len(taus)))*np.nan
    for tind in np.arange(0, max_t):
        mol0s = []
        for tau in taus:
            dr, mol0 = get_dr(rs, PolInd, Nmon, thresh,
                              int(np.floor(tind+tau/2)), Ls)
            mol0s.append(mol0)
        mol0s = np.array(mol0s)
        print(mol0s)
        if not np.isnan(mol0s[0]):
            dRts[tind, :, :] = get_dRt2(rs, mol0s, tind, taus, inds)
    return dRts


def get_dRt0s(rs, Ind, taus, inds):
    Nt = np.shape(rs)[2]
    max_t = Nt - np.max(taus)-1
    dRts = np.zeros((max_t, len(inds), len(taus)))*np.nan
    for tind in np.arange(0, max_t):
        # dr, mol0 = get_dr(rs, PolInd, Nmon, thresh, int(np.floor(tind+tau/2)))
        mol0 = Ind
        if not np.isnan(mol0):
            dRts[tind, :, :] = get_dRt(rs, mol0, tind, taus, inds)
    return dRts


def pol_vs(polatomIDs, vs, rs, maxmon):
    v_pol = []
    for polatomID in polatomIDs:
        dr = get_drPolInd(maxmon, polatomID, rs)
        v_pol.append(
            np.abs(np.dot(vs[int(polatomID-1), :], dr)/np.linalg.norm(dr)))
    return v_pol


def get_drPolInd(maxmon, polatomID, rs):
    dR = np.sum((rs[int(polatomID-1), :] - rs[:maxmon, :])**2, axis=1)
    indmon = np.argmin(dR)
    dr = np.zeros(3)
    if dR[indmon+1] > dR[indmon-1]:
        dr = rs[indmon, :] - rs[indmon-1, :]
    else:
        dr = rs[indmon+1, :] - rs[indmon, :]
    return dr


def get_ribDist(Nrib, seedval, polysome_distfile):
    Dist = np.round(np.loadtxt(polysome_distfile))
    CumDist = np.cumsum(Dist/np.sum(Dist))

    np.random.seed(seedval)
    Nper = np.arange(1, 8)
    Nmade_rib = 0
    Nlist = np.zeros(len(Nper))
    while Nmade_rib < 0.63*Nrib:
        N_i = int(Nper[np.where(np.random.rand() < CumDist)[0][0]])
        if N_i > 1:
            Nlist[N_i] += 1
            Nmade_rib += N_i
    Nfree = Nrib - Nmade_rib
    return Nfree, Nlist


def get_POLforces(filename, Natoms):
    Fs = []
    t = []
    switch = 0
    fs = np.zeros((Natoms, 3))
    f0 = open(filename)
    for line in f0:
        if line.strip()[0:4] == 'Step':
            t.append(float(line.strip().split(' ')[1]))
            Fs = emptydstack(Fs, fs)
            fs = np.zeros((Natoms, 3))
        if switch == 1:
            x = line.strip().split('\t')
            i = int(x[0])-1
            j = int(x[1])-1
            f = np.array([float(x[2]), float(x[3]), float(x[4])])
            fs[i, :] += f
            fs[j, :] -= f
            switch = 0
        if line.strip()[0:11] == 'Pair Forces':
            switch = 1
    f0.close()
    return Fs[:, :, 1:], np.array(t[1:])


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
    r_mag = np.linalg.norm(r)
    costheta = r[2]/r_mag
    theta = np.arccos(costheta)
    sintheta = np.sin(np.arccos(costheta))
    cosphi = r[0]/(r_mag*sintheta)
    sinphi = r[1]/(r_mag*sintheta)
    # Determine quadrant of angle based on signs of cos and sin and
    #    place angle accordingly and convert to degrees
    if cosphi > 0:
        phi = np.arcsin(sinphi)
    elif sinphi < 0:
        phi = -np.arccos(cosphi)
    else:
        phi = np.arccos(cosphi)
    return theta, phi


def write_xyzfile(filename, rs, types, L):  # mols, r_mol, atom_types, time):
    names = ["X", "Ac", "Ag", "Al", "Am", "Ar", "As", "At", "Au", "B", "Ba", "Be", "Bh",
             "Bi", "Bk", "Br", "Ca", "Cd", "Ce", "Cf", "Cl", "Cm", "Co", "Cr", "Cs", "Cu", "Db",
             "Ds", "Dy", "Er", "Es", "Eu", "F", "Fe", "Fm", "Fr", "Ga", "Gd", "Ge", "He", "Hf",
             "Hg", "Ho", "Hs", "I", "In", "Ir", "K", "Kr", "La", "Li", "Lr", "Lu", "Md", "Mg", "Mn",
             "Mo", "Mt", "Na", "Nb", "Nd", "Ne", "Ni", "No", "Np", "Os", "Pa", "Pb",
             "Pd", "Pm", "Po", "Pr", "Pt", "Pu", "Ra", "Rb", "Re", "Rf", "Rg", "Rh", "Rn", "Ru",
             "Sb", "Sc", "Se", "Sg", "Si", "Sm", "Sn", "Sr", "Ta", "Tb", "Tc", "Te", "Th", "Ti", "Tl",
             "Tm", "U", "V", "W", "Xe", "Y", "Yb", "Zn", "Zr"]
    Ntype = len(np.unique(types))
    file_handle = open(filename, 'w')
    N = np.shape(rs)[0]
    Nt = np.shape(rs)[2]
    print(N)
    for t in np.arange(0, Nt):
        file_handle.write(str(N))
        file_handle.write('\nAtoms. Timestep: '+str(0))

        for i in np.arange(0, N):
            file_handle.write('\n'+names[int(types[i])])
            for k in np.arange(0, 3):
                file_handle.write('\t'+str(rs[i, k, t]))
        '''
		for i in np.arange(0,N):
			file_handle.write('\n'+names[int(types[i]+Ntype)])
			for k in np.arange(0,3):
				file_handle.write('\t'+str(rs[i,k,t]-L*np.round(rs[i,k,t]/L)))
		'''
        file_handle.write('\n')
    file_handle.close()


def make_tclfile(tclname, diameters, L):
    f = open(tclname, 'w')
    ntype = len(diameters)
    colors = ["blue", "red", "gray", "orange", "yellow", "tan", "silver", "green", "white",
              "pink", "cyan", "purple", "lime", "mauve", "ochre", "iceblue", "black", "yellow2",
              "yellow3", "green2", "green3", "cyan2", "cyan3", "blue2", "blue3", "violet", "violet2",
              "magenta", "magenta2", "red2", "red3", "orange2", "orange3"]
    names = ["X", "Ac", "Ag", "Al", "Am", "Ar", "As", "At", "Au", "B", "Ba", "Be", "Bh",
             "Bi", "Bk", "Br", "Ca", "Cd", "Ce", "Cf", "Cl", "Cm", "Co", "Cr", "Cs", "Cu", "Db",
             "Ds", "Dy", "Er", "Es", "Eu", "F", "Fe", "Fm", "Fr", "Ga", "Gd", "Ge", "He", "Hf",
             "Hg", "Ho", "Hs", "I", "In", "Ir", "K", "Kr", "La", "Li", "Lr", "Lu", "Md", "Mg", "Mn",
             "Mo", "Mt", "Na", "Nb", "Nd", "Ne", "Ni", "No", "Np", "Os", "Pa", "Pb",
             "Pd", "Pm", "Po", "Pr", "Pt", "Pu", "Ra", "Rb", "Re", "Rf", "Rg", "Rh", "Rn", "Ru",
             "Sb", "Sc", "Se", "Sg", "Si", "Sm", "Sn", "Sr", "Ta", "Tb", "Tc", "Te", "Th", "Ti", "Tl",
             "Tm", "U", "V", "W", "Xe", "Y", "Yb", "Zn", "Zr"]
    i_c = 0
    i_n = 0
    for i_n in np.arange(0, ntype):
        f.write('color Element '+names[i_n]+' '+colors[i_c]+'\n')
        f.write('set natoms [atomselect 0 \"name '+names[i_n]+'\";];\n')
        f.write('$natoms set radius '+str(diameters[i_n]/2)+'\n\n')
        i_n += 1
        i_c += 1
        if i_c > len(colors)-1:
            i_c = 0
    f.write('set cell [pbc set {'+str(L)+' '+str(L)+' '+str(L)+' } -all];\n')
    f.write('pdb box -toggle -center origin -color red;')
    f.close()


ProcPol(args.xyzfile, args.logfile, args.Ndt, args.dt, args.P0, args.f0, args.phi, args.seed, args.Nsteps, args.T, args.epsilon, args.molfile, args.theta,
        args.lef, args.bondangles, args.sigCrowd, args.k_att, args.sig_p, args.N, args.BCs, args.Nmon)
