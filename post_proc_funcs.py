import numpy as np

def make_tclfile(tclname, Nfile, L):
    f=open(tclname,'w')
    Atomfile=np.loadtxt(Nfile)
    ntype=len(Atomfile[:,0])
    diameters=Atomfile[:,2]
    colors=["blue","red","gray","orange","yellow","tan","silver","green","white",
            "pink","cyan","purple","lime","mauve","ochre","iceblue","black","yellow2",
            "yellow3","green2","green3","cyan2","cyan3","blue2","blue3","violet","violet2",
            "magenta","magenta2","red2","red3","orange2","orange3"]
    names=["X","Ac","Ag","Al","Am","Ar","As","At","Au","B","Ba","Be","Bh",
                    "Bi","Bk","Br","Ca","Cd","Ce","Cf","Cl","Cm","Co","Cr","Cs","Cu","Db",
                    "Ds","Dy","Er","Es","Eu","F","Fe","Fm","Fr","Ga","Gd","Ge","He","Hf",
                    "Hg","Ho","Hs","I","In","Ir","K","Kr","La","Li","Lr","Lu","Md","Mg","Mn",
                    "Mo","Mt","Na","Nb","Nd","Ne","Ni","No","Np","Os","Pa","Pb",
                    "Pd","Pm","Po","Pr","Pt","Pu","Ra","Rb","Re","Rf","Rg","Rh","Rn","Ru",
                    "Sb","Sc","Se","Sg","Si","Sm","Sn","Sr","Ta","Tb","Tc","Te","Th","Ti","Tl",
                    "Tm","U","V","W","Xe","Y","Yb","Zn","Zr"]
    i_c=0
    i_n=0
    for i_n in np.arange(0,ntype):
        f.write('color Element '+names[i_n]+' '+colors[i_c]+'\n')
        f.write('set natoms [atomselect 0 \"name '+names[i_n]+'\";];\n')
        f.write('$natoms set radius '+str(diameters[i_n]/2)+'\n\n')
        i_n+=1
        i_c+=1
        if i_c>len(colors)-1:
            i_c=0
    f.write('set cell [pbc set {'+str(L)+' '+str(L)+' '+str(L)+' } -all];\n')
    f.write('pdb box -toggle -center origin -color red;')
    f.close()


def edit_xyzfile(xyzfile):
    names=["X","Ac","Ag","Al","Am","Ar","As","At","Au","B","Ba","Be","Bh",
                "Bi","Bk","Br","Ca","Cd","Ce","Cf","Cl","Cm","Co","Cr","Cs","Cu","Db",
                "Ds","Dy","Er","Es","Eu","F","Fe","Fm","Fr","Ga","Gd","Ge","He","Hf",
                "Hg","Ho","Hs","I","In","Ir","K","Kr","La","Li","Lr","Lu","Md","Mg","Mn",
                "Mo","Mt","Na","Nb","Nd","Ne","Ni","No","Np","Os","Pa","Pb",
                "Pd","Pm","Po","Pr","Pt","Pu","Ra","Rb","Re","Rf","Rg","Rh","Rn","Ru",
                "Sb","Sc","Se","Sg","Si","Sm","Sn","Sr","Ta","Tb","Tc","Te","Th","Ti","Tl",
                "Tm","U","V","W","Xe","Y","Yb","Zn","Zr"]
    f=open(xyzfile,'r')
    f1=open(xyzfile[0:-4]+'_1.xyz','w')
    for x in f:
        line=x.split(' ')
        if len(line)>3:
            newline=names[int(line[0])-1]
            for i in np.arange(1,len(line)):
                newline+='\t'+line[i]
            f1.write(newline)
        else:
            f1.write(x)
    f.close()
    f1.close()

def emptydstack(x,a):
    if len(x)>0:
        x=np.dstack((x,a))
    else:
        x=a
    return x

def process_lmp(lmpfile,dt):
    t=[]
    rs=[]
    vs=[]
    switch=0
    Nval=0
    cnt=0
    tval=0
    l_switch=0
    f=open(lmpfile)
    for x in f:
        print(x.strip())
        if x.strip()=='ITEM: TIMESTEP':
            #t.append(float(f.readline().strip()))
            #print(t[-1])
            tval=1
            pos_data=0
            if switch>0:# and cnt%10==0:
                rs=emptydstack(rs,r)
                vs=emptydstack(vs,v)
                if switch==1:
                    switch=2
        elif tval==1:
            tval=0
            t.append(float(x.strip()))
        elif x.strip()=='ITEM: NUMBER OF ATOMS' and switch==0:
            Nval=1
        elif Nval==1:
            Nval=0
            N=int(x.strip())
        elif x[0:16]=='ITEM: BOX BOUNDS' and switch==0:
            l_switch=1
        elif l_switch==1:
            l_switch=2
            l_line=x.strip().split(' ')
            #print(l_line)
            Lx=-float(l_line[0])+float(l_line[1])
            #print(float(l_line[0]))
  
            #print(Lx)
        elif l_switch==2:
            l_switch=3
            l_line=x.strip().split(' ')
            #print(l_line)
            Ly=-float(l_line[0])+float(l_line[1])
        elif l_switch==3:
            l_switch=0
            l_line=x.strip().split(' ')
            #print(l_line)
            Lz=-float(l_line[0])+float(l_line[1])
            Ls=np.array([Lx,Ly,Lz])
        elif x[0:11]=='ITEM: ATOMS':
            cnt+=1
            pos_data=1
            if switch==0:
                switch=1
                #print('initialize')
                #rs=np.zeros((N,3,1000))
                #vs=np.zeros((N,3,1000))
                #print('initialization complete')
                r=np.zeros((N,3))
                v=np.zeros((N,3))
                types=np.zeros((N))
        elif pos_data==1:
            print(x)
            line=x.strip().split(' ')
            atomid=int(line[0])-1
            x=float(line[2])
            y=float(line[3])
            z=float(line[4])
            #print(line)
            try:
                float(line[5])
            except:
                print(line)
            vx=float(line[5])
            vy=float(line[6])
            vz=float(line[7])
            r[atomid,:]=np.array([x,y,z])
            v[atomid,:]=np.array([vx,vy,vz])
            if switch==1:
                atomtype=int(line[1])-1
                types[atomid]=atomtype
    rs=emptydstack(rs,r)
    vs=emptydstack(vs,v)
    t=np.array(t)*float(dt)
    Natoms=int(N)
    return t, rs, vs, types, Ls, Natoms

def process_lmp2(lmpfile,dt):
    t=[]
    rs=[]
    vs=[]
    switch=0
    Nval=0
    cnt=0
    tval=0
    l_switch=0
    f=open(lmpfile)
    for x in f:
        if x.strip()=='ITEM: TIMESTEP':
            #t.append(float(f.readline().strip()))
            #print(t[-1])
            tval=1
            pos_data=0
            if switch>0: #and cnt%10==0:
                rs=emptydstack(rs,r)
                vs=emptydstack(vs,v)
                if switch==1:
                    switch=2
        elif tval==1:
            tval=0
            t.append(float(x.strip()))
        elif x.strip()=='ITEM: NUMBER OF ATOMS' and switch==0:
            Nval=1
        elif Nval==1:
            Nval=0
            N=int(x.strip())
        elif x[0:16]=='ITEM: BOX BOUNDS' and switch==0:
            l_switch=1
        elif l_switch==1:
            l_switch=2
            l_line=x.strip().split(' ')
            #print(l_line)
            Lx=-float(l_line[0])+float(l_line[1])
            #print(float(l_line[0]))
  
            #print(Lx)
        elif l_switch==2:
            l_switch=3
            l_line=x.strip().split(' ')
            #print(l_line)
            Ly=-float(l_line[0])+float(l_line[1])
        elif l_switch==3:
            l_switch=0
            l_line=x.strip().split(' ')
            #print(l_line)
            Lz=-float(l_line[0])+float(l_line[1])
            Ls=np.array([Lx,Ly,Lz])
        elif x[0:11]=='ITEM: ATOMS':
            cnt+=1
            pos_data=1
            if switch==0:
                switch=1
                #print('initialize')
                #rs=np.zeros((N,3,1000))
                #vs=np.zeros((N,3,1000))
                #print('initialization complete')
                r=np.zeros((N,3))
                v=np.zeros((N,3))
                types=np.zeros((N))
        elif pos_data==1:
            line=x.strip().split(' ')
            atomid=int(line[0])-1
            x=float(line[2])
            y=float(line[3])
            z=float(line[4])
            r[atomid,:]=np.array([x,y,z])
            if switch==1:
                atomtype=int(line[1])-1
                types[atomid]=atomtype
    rs=emptydstack(rs,r)
    t=np.array(t)*float(dt)
    Natoms=int(N)
    return t, rs, types, Ls, Natoms

def subtract_com(arr,types, Nfile):
    Atomfile=np.loadtxt(Nfile)
    masses=Atomfile[:,1]
    type_arr=[]
    for typeA in types:
        type_arr.append(int(typeA))
    type_arr=np.array(type_arr)
    massval=masses[types]
    r_com=np.sum(arr*massval[:,None,None],axis=0)/np.sum(massval)
    arr=arr-r_com[None,:,:]
    return arr

def subtract_com2(arr,types, masses):
    type_arr=[]
    for typeA in types:
        type_arr.append(int(typeA))
    type_arr=np.array(type_arr)
    massval=masses[types]
    r_com=np.sum(arr*massval[:,None,None],axis=0)/np.sum(massval)
    arr=arr-r_com[None,:,:]
    return arr

def subtract_comSticky(arr,types, masses, attractions):
    type_arr=[]
    for typeA in types:
        type_arr.append(int(typeA))
    type_arr=np.array(type_arr)
    massval=masses[types]*attractions
    r_com=np.sum(arr*massval[:,None,None],axis=0)/np.sum(massval)
    arr=arr-r_com[None,:,:]
    return arr


def make_MSDMFD(t, rs, Natoms, Ndt):
    inds=np.unique(np.round(np.logspace(np.log10(1),np.log10(len(t)-2),Ndt)))
    dts=inds*float(t[1]-t[0])
    MSD=np.zeros((Natoms,len(inds)))
    MFD=np.zeros((Natoms,len(inds)))
    for i in np.arange(0,len(inds)):
        ind=int(inds[i])
        drs=rs[:,:,0:-ind]-rs[:,:,ind:]
        MSD[:,i]=np.mean(np.sum(drs**2,axis=1),axis=1)
        MFD[:,i]=np.mean(np.sum(drs**4,axis=1),axis=1)
    return dts, MSD, MFD

def make_rs2(rs, n_steps):
    rs2 = np.empty((np.shape(rs)[0], np.shape(rs)[1], len(n_steps), int(np.ceil(np.shape(rs)[2]/len(n_steps)))))*np.nan
    for i in np.arange(0,np.ceil(np.shape(rs)[2]/len(n_steps))):
        rs2[:,:,:,i] = rs[:,:,i*len(n_steps):(i+1)*len(n_steps)]
    return rs2

def make_MSDMFD2(t, rs):
    '''
    dR = rs[:,:,1:] - rs[:,:,0:-1]
    dt = np.round(t[1:] - t[0:-1], decimals=6)
    dt_uni = np.unique(dt)
    dR2 = np.ones((np.shape(dR)[0],np.shape(dR)[1], len(dt_uni), int(np.ceil(np.shape(dR)[2]/len(dt_uni))+2)))*np.nan
    dt2 = np.ones((np.shape(dR)[0],np.shape(dR)[1], len(dt_uni), int(np.ceil(np.shape(dR)[2]/len(dt_uni))+2)))*np.nan
    ind_arr = np.zeros(len())
    ind =0
    for i in np.arange(0,np.shape(dR)[2]):
        dR2[:,:,np.where(dt[i]==dt_uni)[0][0],ind] = dR[:,:,i]
        if dt[i]==dt_uni[0]:
            ind+=1
    drs = dR2
    MSD = np.nanmean(np.nansum(dR2**2,axis=1),axis=2)
    MFD = np.nanmean(np.nansum(dR2**4,axis=1),axis=2)
    '''
    dR = rs[:,:,1:] - rs[:,:,0:-1]
    dt = np.round(t[1:] - t[0:-1], decimals=6)
    dt_uni = np.unique(dt)
    MSD = np.zeros((np.shape(dR)[0], len(dt_uni)))
    MFD = np.zeros((np.shape(dR)[0], len(dt_uni)))
    N = np.zeros(len(dt_uni))
    for i in np.arange(0,np.shape(dR)[2]):
        print(np.where(dt[i]==dt_uni)[0][0])
        MSD[:,np.where(dt[i]==dt_uni)[0][0]]+= np.nansum(dR[:,:,i]**2,axis=1)
        MFD[:,np.where(dt[i]==dt_uni)[0][0]]+= np.nansum(dR[:,:,i]**4,axis=1)
        N[np.where(dt[i]==dt_uni)[0][0]]+=1
    MSD/=N
    MFD/=N

    return dt_uni, MSD, MFD

def make_VAC2(t, vs):
    '''
    vdv = vs[:,:,1:] * vs[:,:,0:-1]
    dt = np.round(t[1:] - t[0:-1], decimals=6)
    dt_uni = np.unique(dt)
    VdV2 = np.ones((np.shape(vdv)[0],np.shape(vdv)[1], len(dt_uni), int(np.ceil(np.shape(vdv)[2]/len(dt_uni))+2)))*np.nan
    ind =0
    for i in np.arange(0,np.shape(vdv)[2]):
        #print(np.where(dt[i]==dt_uni)[0][0])
        VdV2[:,:,np.where(dt[i]==dt_uni)[0][0],ind] = vdv[:,:,i]
        if dt[i]==dt_uni[0]:
            ind+=1
    VAC = np.nanmean(np.nansum(VdV2,axis=1),axis=2)
    '''
    vdv = vs[:,:,1:] * vs[:,:,0:-1]
    dt = np.round(t[1:] - t[0:-1], decimals=6)
    dt_uni = np.unique(dt)
    VAC = np.zeros((np.shape(vdv)[0],len(dt_uni)))
    N = np.zeros(len(dt_uni))
    for i in np.arange(0,np.shape(vdv)[2]):
        #print(np.where(dt[i]==dt_uni)[0][0])
        VAC[:,np.where(dt[i]==dt_uni)[0][0]] += np.nansum(vdv[:,:,i],axis=1)
        N[np.where(dt[i]==dt_uni)[0][0]] +=1
    VAC/=N
    return dt_uni, VAC

def make_ISF2(t, rs, ks, Natoms):
    '''
    dt = np.round(t[1:] - t[0:-1], decimals=6)
    dt_uni = np.unique(dt)
    ISF=np.zeros((Natoms,len(dt_uni),len(ks)))
    dR = rs[:,:,1:] - rs[:,:,0:-1]
    dR2 = np.ones((np.shape(dR)[0],np.shape(dR)[1], len(dt_uni), int(np.ceil(np.shape(dR)[2]/len(dt_uni))+2)))*np.nan
    dt2 = np.ones((np.shape(dR)[0],np.shape(dR)[1], len(dt_uni), int(np.ceil(np.shape(dR)[2]/len(dt_uni))+2)))*np.nan
    ind =0
    for i in np.arange(0,np.shape(dR)[2]):
        #print(np.where(dt[i]==dt_uni)[0][0])
        dR2[:,:,np.where(dt[i]==dt_uni)[0][0],ind] = dR[:,:,i]
        if dt[i]==dt_uni[0]:
            ind+=1
    drs = dR2
    for j in np.arange(0,len(ks)):
        k=ks[j]/np.sqrt(len(rs[0,:,0])) # Given isotropy take vector 1/sqrt(3)(i, j, k)
        ISF[:,:,j]=np.nanmean(np.exp(j*k*np.nansum(drs,axis=1)).real,axis=2)
    '''
    dt = np.round(t[1:] - t[0:-1], decimals=6)
    dt_uni = np.unique(dt)
    ISF=np.zeros((Natoms,len(dt_uni),len(ks)))
    dR = rs[:,:,1:] - rs[:,:,0:-1]
    N = np.zeros(len(dt_uni))
    for i in np.arange(0,np.shape(dR)[2]):
        for j in np.arange(0,len(ks)):
            k=ks[j]/np.sqrt(len(rs[0,:,0]))
            ISF[:,np.where(dt[i]==dt_uni)[0][0],j] += np.exp(1j*k*np.nansum(dR[:,:,i])).real
            N[np.where(dt[i]==dt_uni)[0][0]] += 1
    for j in np.arange(0,len(ks)):
        ISF[:,:,j]/=N
    return dt_uni, ISF

def make_VAC(t, vs, Natoms, Ndt):
	inds=np.unique(np.round(np.logspace(np.log10(1),np.log10(len(t)-2),Ndt)))
	dts=inds*float(t[1]-t[0])
	VAC=np.zeros((Natoms,len(inds)))
	for i in np.arange(0,len(inds)):
		ind=int(inds[i])
		vdv=vs[:,:,0:-ind]*vs[:,:,ind:]
		VAC[:,i]=np.mean(np.sum(vdv,axis=1),axis=1)
	return dts, VAC

def make_ISF(t, rs, ks, Natoms, Ndt):
	inds=np.unique(np.round(np.logspace(np.log10(1),np.log10(len(t)-2),Ndt)))
	dts=inds*float(t[1]-t[0])
	ISF=np.zeros((Natoms,len(dts),len(ks)))
	for i in np.arange(0,len(inds)):
		ind=int(inds[i])
		drs=rs[:,:,0:-ind]-rs[:,:,ind:]
		for j in np.arange(0,len(ks)):
			k=ks[j]/np.sqrt(len(rs[0,:,0]))
			ISF[:,i,j]=np.mean(np.exp(1j*k*np.sum(drs,axis=1)).real,axis=1)
	return dts, ISF

def get_taus(dts, ISF):
    Ntype=len(ISF[:,0,0])
    Nk=len(ISF[0,0,:])
    taus=np.zeros((Ntype,Nk))
    tau_inds=np.zeros((Ntype,Nk))
    for i in np.arange(0,Ntype):
        for j in np.arange(0,Nk):
            inds=np.where(ISF[i,:,j]<np.exp(-1))
            if np.shape(inds)[0]>0 and np.shape(inds)[1]>0:
                tau_inds[i,j]=int(np.where(ISF[i,:,j]<np.exp(-1))[0][0])
                taus[i,j]=dts[int(tau_inds[i,j])]
            else:
                tau_inds[i,j]=np.nan
                taus[i,j]=np.nan

    return taus, tau_inds

def make_Pr(dts,rs,tau_inds,bin_edges, types):
    # bin_edges has size Ntau, Nbins+1
    tau_inds=tau_inds.astype(int)
    ntypes=np.max(types)+1
    Ntau=len(tau_inds[0,:])
    Pr=np.ones((ntypes,Ntau,len(bin_edges[0,:])-1))*np.nan
    for j in np.arange(0,Ntau):
    	for i in np.arange(0,np.max(types)+1):
    		ind=tau_inds[i,j]
    		if (not np.isnan(ind))and ind>0:
    			bin_arr=bin_edges[j,:]
    			type_ind=np.where(types==i)[0]
    			#print(ind)
    			#print(np.shape(rs[type_ind,:,0:-ind]))
    			#print(np.shape(rs[type_ind,:,0:-ind].flatten()))
    			#print(np.shape(rs[type_ind,:,ind]))
    			drs=np.sqrt(np.sum((rs[type_ind,:,0:-ind]-rs[type_ind,:,ind:])**2,axis=1)).flatten()
    			Pr[i,j,:], bin_es=np.histogram(drs,bins=bin_arr,density=True)
    return Pr

def make_DAC(t, Ndt, rs, tau_inds, types):
	inds=np.unique(np.round(np.logspace(np.log10(1),np.log10(len(t)-2),Ndt)))
	dts=inds*float(t[1]-t[0])
	ntypes=np.max(types)+1
	Ntau=len(tau_inds[0,:])
	DAC=np.empty((ntypes,len(inds),Ntau))*np.nan
	for k in np.arange(0,ntypes):
		for j in np.arange(0,Ntau):	
			tau_ind=tau_inds[k,j]
			if (not np.isnan(tau_ind)) and tau_ind>0:		
				tau_ind=int(tau_ind)
				for i in np.arange(0,len(inds)-tau_ind):
					ind=int(inds[i])
					type_inds=np.where(types==k)[0]
					drs1=rs[type_inds,:,tau_ind:-ind]-rs[type_inds,:,:-(tau_ind+ind)]
					drs2=rs[type_inds,:,tau_ind+ind:]-rs[type_inds,:,ind:-tau_ind]
					DAC[k,i,j]=np.mean(np.mean(np.sum(drs1*drs2,axis=1),axis=1))
	return DAC

def post_proc(lmpfile, datfile, Nfile, Ndt, dt):
	t, rs, vs, types, Ls, Natoms=process_lmp(lmpfile,dt)
	Atomfile=np.loadtxt(Nfile)
	diameters=Atomfile[:,2]
	## Select ks based of of relevant length scales in system
	l_0s=np.array([np.min(diameters), np.max(diameters),1,float(Ls[0])])
	ks=2*np.pi/l_0s
	types=types.astypes(int)
	ntypes=np.max(types)+1
	Nbins=int(40.0)
	bin_edges=np.zeros((len(ks),Nbins+1))
	for i in np.arange(0,len(l_0s)):
		bin_edges[i,:]=np.arange(0,2.0*l_0s[i],2.0*l_0s[i]/float(Nbins+1))
	bin_centers=(bin_edges[:,1:]+bin_edges[:,:-1])/2.0
	## Subtract off center of mass of velocities and positons
	rs=subtract_com(rs,types, Nfile)
	vs=subtract_com(vs,types, Nfile)
	## Calculate relevant statistics
	dts, MSD, MFD=make_MSDMFD(t, rs, Natoms, Ndt)
	dts, VAC=make_VAC(t, vs, Natoms, Ndt)
	dts, ISF=make_ISF(t, rs, ks, Natoms, Ndt)
	## Average over atom type
	VACo=np.zeros((ntypes,len(dts)))
	MSDo=np.zeros((ntypes,len(dts)))
	MFDo=np.zeros((ntypes,len(dts)))
	ISFo=np.zeros((ntypes,len(dts),len(ks)))
	for i in np.arange(0,ntypes):
		VACo[i,:]=np.mean(VAC[np.where(types==i)[0],:],axis=0)
		MSDo[i,:]=np.mean(MSD[np.where(types==i)[0],:],axis=0)
		MFDo[i,:]=np.mean(MFD[np.where(types==i)[0],:],axis=0)
		ISFo[i,:,:]=np.mean(ISF[np.where(types==i)[0],:,:],axis=0)
	## Calculate relaxation timescales
	taus, tau_inds = get_taus(dts, ISFo)
	## Calculate P(dr(tau))
	Pr=make_Pr(dts,rs,tau_inds,bin_edges, types)
	## Calculate DAC
	DAC=make_DAC(t, Ndt, rs, tau_inds, types)
	## Calculate non-gaussian parameter
	A2o=3.0*MFDo/(5.0*MSDo**2.0)-1.0
	## Output statistics into datfile
	f=open(datfile,'w')
	for i in np.arange(0,len(dts)):
		f.write(str(dts[i]))
		for j in np.arange(0,ntypes):
			f.write('\t'+str(MSDo[j,i]))
		for j in np.arange(0,ntypes):
			f.write('\t'+str(MFDo[j,i]))
		for j in np.arange(0,ntypes):
			f.write('\t'+str(A2o[j,i]))
		for j in np.arange(0,ntypes):
			f.write('\t'+str(VACo[j,i]))
		f.write('\n')
	f.close()

	## Write output for isf file
	f=open(datfile[:-4]+'.isf','w')
	type_arr=np.arange(0,ntypes)
	f.write(str(np.nan))
	for i in np.arange(0,len(l_0s)):
		for j in np.arange(0,ntypes):
			f.write('\t'+str(j))
	f.write('\n')
	f.write(str(0))
	for i in np.arange(0,len(l_0s)):
		for j in np.arange(0,ntypes):
			f.write('\t'+str(l_0s[i]))
	f.write('\n')
	f.write(str(0))
	for i in np.arange(0,len(l_0s)):
		for j in np.arange(0,ntypes):
			f.write('\t'+str(ks[i]))
	f.write('\n')
	f.write(str(0))
	for i in np.arange(0,len(l_0s)):
		for j in np.arange(0,ntypes):
			f.write('\t'+str(taus[j,i]))
	f.write('\n')
	for k in np.arange(0,len(ISFo[0,:,0])):
		f.write(str(dts[k]))
		for i in np.arange(0,len(l_0s)):
			for j in np.arange(0,ntypes):
				f.write('\t'+str(ISFo[j,k,i]))
		f.write('\n')
	f.close()

	## Write output for dac file
	f=open(datfile[:-4]+'.dac','w')
	type_arr=np.arange(0,ntypes)
	f.write(str(np.nan))
	for i in np.arange(0,len(l_0s)):
		for j in np.arange(0,ntypes):
			f.write('\t'+str(j))
	f.write('\n')
	f.write(str(0))
	for i in np.arange(0,len(l_0s)):
		for j in np.arange(0,ntypes):
			f.write('\t'+str(l_0s[i]))
	f.write('\n')
	f.write(str(0))
	for i in np.arange(0,len(l_0s)):
		for j in np.arange(0,ntypes):
			f.write('\t'+str(ks[i]))
	f.write('\n')
	f.write(str(0))
	for i in np.arange(0,len(l_0s)):
		for j in np.arange(0,ntypes):
			f.write('\t'+str(taus[j,i]))
	f.write('\n')
	for k in np.arange(0,len(DAC[0,:,0])):
		f.write(str(dts[k]))
		for i in np.arange(0,len(l_0s)):
			for j in np.arange(0,ntypes):
				f.write('\t'+str(DAC[j,k,i]))
		f.write('\n')
	f.close()

	## Write output for pdr file
	f=open(datfile[:-4]+'.pdr','w')
	type_arr=np.arange(0,ntypes)
	f.write(str(np.nan))
	for i in np.arange(0,len(l_0s)):
		for j in np.arange(0,ntypes):
			f.write('\t'+str(j))
		f.write('\t'+str(np.nan))
	f.write('\n')
	f.write(str(np.nan))
	for i in np.arange(0,len(l_0s)):
		for j in np.arange(0,ntypes):
			f.write('\t'+str(l_0s[i]))
		f.write('\t'+str(np.nan))
	f.write('\n')
	f.write(str(np.nan))
	for i in np.arange(0,len(l_0s)):
		for j in np.arange(0,ntypes):
			f.write('\t'+str(ks[i]))
		f.write('\t'+str(np.nan))
	f.write('\n')
	f.write(str(np.nan))
	for i in np.arange(0,len(l_0s)):
		for j in np.arange(0,ntypes):
			f.write('\t'+str(taus[j,i]))
		f.write('\t'+str(np.nan))
	f.write('\n')
	for k in np.arange(0,len(Pr[0,0,:])):
		f.write(str(dts[k]))
		for i in np.arange(0,len(l_0s)):
			if i==0:
				f.write(str(bin_centers[i,k]))
			else:
				f.write('\t'+str(bin_centers[i,k]))
			for j in np.arange(0,ntypes):
				f.write('\t'+str(Pr[j,i,k]))
		f.write('\n')
	f.close()



def make_molfile(molfile, rs, atomtypes, blist, blengths, balist, thetas, dalist, phis):
    '''
    Function that creates a lammps file.mol for a molecule in lammps (molecular sim)
    
    Inputs:
        molfile         name of file to be created
        rs              array of position vectors [shape: Nx#dim] (array of floats)
        atomtypes       vector of type of each atom corresponding with rs (array of ints)
        blist           numpy array of bond pair [shape: Nbondx2] (array of ints)
        blengths        vector of bond lengths corresponding to blist  (vector of floats)
        balist          numpy array of bond angle sets [shape: Nbax3] (array of ints)
        thetas          vector of bond angles corresponding to balist  (vector of floats)
        dalist          numpy array of dihedral angle sets [shape: Ndax4] (array of ints)
        phis            vector of dihedral angles corresponding to dalist (vector of floats)
    '''
    f=open(molfile,'w')
    Natoms=int(np.shape(rs)[0])
    
    
    f.write('# Chiral Rod\n\n')
    
    f.write(str(Natoms)+' atoms\n')
    if np.shape(blist)[0]>0:
        f.write(str(np.shape(blist)[0])+' bonds \n')
    if np.shape(balist)[0]>0:
        f.write(str(np.shape(balist)[0])+' angles \n')
    if np.shape(dalist)[0]>0:
        f.write(str(np.shape(dalist)[0])+' dihedrals \n')
    f.write('\n')
    
    f.write('Coords\n\n')
    for i in np.arange(0,Natoms):
        f.write(str(i+1)+' '+str(rs[i,0])+' '+str(rs[i,1])+' '+str(rs[i,2])+'\n')
    f.write('\n')
    
    f.write('Types\n\n')
    for i in np.arange(0,Natoms):
        f.write(str(i+1)+' '+str(int(atomtypes[i]+1))+'\n')
    f.write('\n')
    
    if np.shape(blist)[0]>0:
        blength_vals=np.unique(blengths)
        f.write('Bonds \n\n')
        for i in np.arange(0,np.shape(blist)[0]):
            pair=blist[i,:]+1
            blength=blengths[i]
            type_val=np.where(blength_vals==blength)[0][0]+1
            f.write(str(i+1)+' '+str(type_val)+' '+str(int(pair[0]))+' '+str(int(pair[1]))+'\n')
        f.write('\n')
    if np.shape(balist)[0]>0:
        f.write('Angles \n\n')
        theta_vals=np.unique(thetas)
        for i in np.arange(0,np.shape(balist)[0]):
            trio=balist[i,:]+1
            theta=thetas[i]
            type_val=np.where(theta_vals==theta)[0][0]+1
            f.write(str(i+1)+' '+str(type_val)+' '+str(int(trio[0]))+' '+str(int(trio[1]))+' '+str(int(trio[2]))+'\n')
        f.write('\n')
    if np.shape(dalist)[0]>0:
        f.write('Dihedrals \n\n')
        phi_vals=np.unique(phis)
        for i in np.arange(0,np.shape(dalist)[0]):
            quad=dalist[i,:]+1
            phi=phis[i]
            type_val=np.where(phi_vals==phi)[0][0]+1
            f.write(str(i+1)+' '+str(type_val)+' '+str(int(quad[0]))+' '+str(int(quad[1]))+' '+str(int(quad[2]))+' '+str(int(quad[3]))+'\n')
        f.write('\n')
    '''    
    if len(blengths)>0:
        # Get Lists for Special Bonds
        list12, list13, list14= get_SpecialLists(blist)
        f.write('Special Bond Counts\n\n')
        for i in np.arange(0,Natoms):
            f.write(str(int(i+1))+' '+str(len(list12[i]))+' '+str(int(len(list13[i])*len(balist)))+' '+str(int(len(list14[i])*len(dalist)))+'\n')
        f.write('\n')
        f.write('Special Bonds\n\n')
        for i in np.arange(0,Natoms):
            bond_str=str(int(i+1))
            for j in list12[i]:
                # Shift indices 1 for going from python to LAMMPS
                j+=1
                bond_str+=' '+str(int(j))
            if len(balist)>0:
                for j in list13[i]:
                    # Shift indices 1 for going from python to LAMMPS
                    j+=1
                    bond_str+=' '+str(int(j))
            if len(dalist)>0:
                for j in list14[i]:
                    # Shift indices 1 for going from python to LAMMPS
                    j+=1
                    bond_str+=' '+str(int(j))
            bond_str+='\n'
            f.write(bond_str)
    '''       
    f.write('\n')
    
    f.close()    
    
def get_SpecialLists(blist):
    Natoms=int(np.max(blist)+1)
    list12=[[] for i in range(Natoms)]
    for pair in blist:
        i=int(pair[0])
        j=int(pair[1])
        list12[i].append(j)
        list12[j].append(i)
    list13=[[] for i in range(Natoms)]
    list14=[[] for i in range(Natoms)]
    i=0
    for i in np.arange(0,Natoms):
        for j in list12[i]:
            next_neighbors=list12[j]
            for k in next_neighbors:
                if k not in list12[i] and k not in list13[i] and k !=i:
                    list13[i].append(k)
    for i in np.arange(0,Natoms):
        for j in list13[i]:
            next_neighbors=list12[j]
            for k in next_neighbors:
                if k not in list12[i] and k not in list13[i] and k not in list14[i] and k !=i:
                    list14[i].append(k)
    return list12, list13, list14