import numpy as np
import h5py as h5

from datetime import datetime, timedelta
from axis_hurr import *
from obs_subs import *
from matrix_reshape import * # including parameter modelpath

stime=datetime(2018,9,11,3,24,20) # real value
etime=datetime(2018,9,11,6)
#bfilename=modelpath+'fwdrun/'+'Modelout-{:%Y%m%d%H}-{:%Y%m%d%H}.h5'\
#    .format(stime, etime)
bfilename=modelpath+'fwdrun/Ideal-35000.h5'

bot=1
top=20+1 # 20 km

# GFS first guess
"""
def get_fg():
    with h5.File(bfilename, 'r') as bfile:
        bu9 = bfile['u9st'][:,:,0]
        bv9 = bfile['v9st'][:,:,0]
        bw9 = bfile['w9st'][:,:,0]
        bt9 = bfile['t9st'][:,:,0]
        bp9 = bfile['p9st'][:,:,0]
        bqv9 = bfile['qv9st'][:,:,0]
        bql9 = bfile['ql9st'][:,:,0]
    return [bu9 , bv9 , bw9 , bt9 , bp9 , bqv9 , bql9]
"""
def get_fg():
    with h5.File(bfilename, 'r') as bfile:
        bu9 = bfile['u9'][:,:]
        bv9 = bfile['v9'][:,:]
        bw9 = bfile['w9'][:,:]
        bt9 = bfile['t9'][:,:]
        bp9 = bfile['p9'][:,:]
        bqv9 = bfile['qv9'][:,:]
        bql9 = bfile['ql9'][:,:]
    return [bu9 , bv9 , bw9 , bt9 , bp9 , bqv9 , bql9]

def get_shapes():
    varlist=get_fg()
    shapes=[]
    for ivar in varlist:
        shapes.append(ivar.shape)
    return shapes
shapes=get_shapes()
    
def get_moderr():
    berrorname=modelpath+'/error-model/ModelErrors.h5'
    #berrorname=modelpath+'/ideal/ModelError.h5'
    with h5.File(berrorname, 'r') as berror:
        uerror=berror['uerror'][()].mean()
        verror=berror['verror'][()].mean()
        werror=berror['werror'][()].mean()
        terror=berror['terror'][()].mean()
        perror=berror['perror'][()].mean()
        qverror=berror['qverror'][()].mean()
        qlerror=berror['qlerror'][()].mean()

    #moderr=np.array([uerror ,verror ,werror ,terror ,perror ,qverror ,qlerror])
    #return moderr
    return [uerror ,verror ,werror ,terror ,perror ,qverror ,qlerror]

def get_precon():
    varerr=get_moderr()
    precon_nd=[]
    for ivar in range(len(shapes)):
        ipre=np.ones(shapes[ivar])
        ipre[1:-1,bot:top]=np.sqrt(varerr[ivar])
        #ipre[1:-1,bot:top]=np.sqrt(varerr[ivar][bot-1:top-1])
        precon_nd.append(ipre)
    buf, precon=nd2oned(precon_nd)
    return precon

def cost_mod(avar1d):
    Jb=0  # the result to be returned

    bu9 , bv9 , bw9 , bt9 , bp9 , bqv9 , bql9 = get_fg()
    au9 , av9 , aw9 , at9 , ap9 , aqv9 , aql9 = oned2nd(shapes, avar1d)

    uerror, verror, werror, terror, perror, qverror, qlerror =get_moderr()

    for ilev in range(bot, top):
        Jb+=((( au9- bu9)[1:-1,ilev]**2)/ uerror).sum() + \
            ((( av9- bv9)[1:-1,ilev]**2)/ verror).sum() + \
            ((( aw9- bw9)[1:-1,ilev]**2)/ werror).sum() + \
            ((( at9- bt9)[1:-1,ilev]**2)/ terror).sum() + \
            ((( ap9- bp9)[1:-1,ilev]**2)/ perror).sum() + \
            (((aqv9-bqv9)[1:-1,ilev]**2)/qverror).sum() + \
            (((aql9-bql9)[1:-1,ilev]**2)/qlerror).sum()

    return Jb*0.5

sstep=10
def cost_obs(avar1d, obs_dict):
    """
    input: 
       avar1d:  1-dimensional analysis variables
       obs_dict: ["inst", "obstime", "fname"]
    """
    Jo=0  # the result to be returned
    au9 , av9 , aw9 , at9 , ap9 , aqv9 , aql9 = oned2nd(shapes, avar1d)
    
    # starting from today
    H=hurr_vars
    hurr_in()
    param_setup()
    hurr_initial()

    double_timer()
    H.u9 = au9.copy()
    H.v9 = av9.copy()
    H.w9 = aw9.copy()
    H.t9 = at9.copy()
    H.p9 = ap9.copy()
    H.qv9 = aqv9.copy()
    H.ql9 = aql9.copy()

    H.u19[()]  = H.u9[ 1:  ,1:-1]
    H.v19[()]  = H.v9[ 1:-1,1:-1]
    H.w19[()]  = H.w9[ 1:-1,1:  ]
    H.t19[()]  = H.t9[ 1:-1,1:-1]
    H.p19[()]  = H.p9[ 1:-1,1:-1]
    H.qv19[()] = H.qv9[1:-1,1:-1]
    H.ql19[()] = H.ql9[1:-1,1:-1]

    estep=int((obs_dict['obstime']-stime).total_seconds()/H.dt)
    for H.itime in range(sstep+1, sstep+estep):
        hurr_fwd()

    if obs_dict['inst']=='AMSR2':
        Jo+=cost_amsr2(H.t9, H.p9, H.qv9, H.ql9, H.pn, obs_dict)
    else:
        Jo+=cost_atms(H.t9.copy(), obs_dict)

    return Jo*0.5

def cost(avar1d):
    precon=get_precon()
    ax=avar1d*precon # convert from L-1*x to x

    obs_dicts=get_obs_dicts()    
    Jobs=0
    for obs_dict in obs_dicts: 
        Jobs+=cost_obs(ax, obs_dict)
    res=Jobs+cost_mod(ax)
    return res

def grad_cost_mod(avar1d):
    bu9 , bv9 , bw9 , bt9 , bp9 , bqv9 , bql9 = get_fg()
    au9 , av9 , aw9 , at9 , ap9 , aqv9 , aql9 = oned2nd(shapes, avar1d)
    uerror, verror, werror, terror, perror, qverror, qlerror = get_moderr()

    du  = np.zeros( bu9.shape, float)
    dv  = np.zeros( bv9.shape, float)
    dw  = np.zeros( bw9.shape, float)
    dt  = np.zeros( bt9.shape, float)
    dp  = np.zeros( bp9.shape, float)
    dqv = np.zeros(bqv9.shape, float)
    dql = np.zeros(bql9.shape, float)
    
    for ilev in range(bot, top):
        du[1:-1,ilev] = (au9-bu9)[1:-1,ilev]/uerror
        dv[1:-1,ilev] = (av9-bv9)[1:-1,ilev]/verror
        dw[1:-1,ilev] = (aw9-bw9)[1:-1,ilev]/werror
        dt[1:-1,ilev] = (at9-bt9)[1:-1,ilev]/terror
        dp[1:-1,ilev] = (ap9-bp9)[1:-1,ilev]/perror
        dqv[1:-1,ilev] = (aqv9-bqv9)[1:-1,ilev]/qverror
        dql[1:-1,ilev] = (aql9-bql9)[1:-1,ilev]/qlerror

    dvarnd=[du , dv , dw , dt , dp , dqv , dql]
    buf, dvar1d=nd2oned(dvarnd)

    return dvar1d
    
def grad_cost_obs(avar1d, obs_dict):
    bu9 , bv9 , bw9 , bt9 , bp9 , bqv9 , bql9 = get_fg()    
    au9 , av9 , aw9 , at9 , ap9 , aqv9 , aql9 = oned2nd(shapes, avar1d)

    H=hurr_vars

    hurr_in()
    param_setup()
    hurr_initial()
    
    double_timer()
    H.u9 = au9.copy()
    H.v9 = av9.copy()
    H.w9 = aw9.copy()
    H.t9 = at9.copy()
    H.p9 = ap9.copy()
    H.qv9 = aqv9.copy()
    H.ql9 = aql9.copy()

    H.u19[()]  = H.u9[ 1:  ,1:-1]
    H.v19[()]  = H.v9[ 1:-1,1:-1]
    H.w19[()]  = H.w9[ 1:-1,1:  ]
    H.t19[()]  = H.t9[ 1:-1,1:-1]
    H.p19[()]  = H.p9[ 1:-1,1:-1]
    H.qv19[()] = H.qv9[1:-1,1:-1]
    H.ql19[()] = H.ql9[1:-1,1:-1]

    estep=int((obs_dict['obstime']-stime).total_seconds()/H.dt)

    # fwd run and storing state variables
    u9st   = np.zeros([au9.shape[0], au9.shape[1], estep], float)
    v9st   = np.zeros([av9.shape[0], av9.shape[1], estep], float)
    w9st   = np.zeros([aw9.shape[0], aw9.shape[1], estep], float)
    t9st   = np.zeros([at9.shape[0], at9.shape[1], estep], float)
    p9st   = np.zeros([ap9.shape[0], ap9.shape[1], estep], float)
    qv9st  = np.zeros([aqv9.shape[0], aqv9.shape[1], estep], float)
    ql9st  = np.zeros([aql9.shape[0], aql9.shape[1], estep], float)

    u19st  = np.zeros([H.u19.shape[0], H.u19.shape[1], estep], float)
    v19st  = np.zeros([H.v19.shape[0], H.v19.shape[1], estep], float)
    w19st  = np.zeros([H.w19.shape[0], H.w19.shape[1], estep], float)
    t19st  = np.zeros([H.t19.shape[0], H.t19.shape[1], estep], float)
    p19st  = np.zeros([H.p19.shape[0], H.p19.shape[1], estep], float)
    qv19st = np.zeros([H.qv19.shape[0], H.qv19.shape[1], estep], float)
    ql19st = np.zeros([H.ql19.shape[0], H.ql19.shape[1], estep], float)

    for H.itime in range(sstep+1, sstep+estep):
        u9st[:,:,H.itime-sstep-1]   = H.u9.copy()
        v9st[:,:,H.itime-sstep-1]   = H.v9.copy()
        w9st[:,:,H.itime-sstep-1]   = H.w9.copy()
        t9st[:,:,H.itime-sstep-1]   = H.t9.copy()
        p9st[:,:,H.itime-sstep-1]   = H.p9.copy()
        qv9st[:,:,H.itime-sstep-1]  = H.qv9.copy()
        ql9st[:,:,H.itime-sstep-1]  = H.ql9.copy()
        u19st[:,:,H.itime-sstep-1]  = H.u19.copy()
        v19st[:,:,H.itime-sstep-1]  = H.v19.copy()
        w19st[:,:,H.itime-sstep-1]  = H.w19.copy()
        t19st[:,:,H.itime-sstep-1]  = H.t19.copy()
        p19st[:,:,H.itime-sstep-1]  = H.p19.copy()
        qv19st[:,:,H.itime-sstep-1] = H.qv19.copy()
        ql19st[:,:,H.itime-sstep-1] = H.ql19.copy()
        
        hurr_fwd()

    H.u[:,:]=0
    H.v[:,:]=0
    H.w[:,:]=0
    H.t[:,:]=0  #dT_ATMS.copy()+dT_AMSR2.copy()
    H.p[:,:]=0  #dP
    H.qv[:,:]=0 #dQV
    H.ql[:,:]=0 #dQL
    H.u1[:,:]=0
    H.v1[:,:]=0
    H.w1[:,:]=0
    H.t1[:,:]=0
    H.p1[:,:]=0
    H.qv1[:,:]=0
    H.ql1[:,:]=0
    
    if obs_dict['inst']=='AMSR2':
        dT, dP, dQV, dQL=grad_cost_amsr2(H.t9, H.p9, H.qv9, H.ql9, H.pn, obs_dict)
        H.t[ :,:]+=dT
        H.p[ :,:]+=dP
        H.qv[:,:]+=dQV
        H.ql[:,:]+=dQL
    else:
        dT=grad_cost_atms(H.t9.copy(), obs_dict)
        H.t[:,:]+=dT

    for itime in range(estep-2, -1, -1): 
        H.u9[:,:] = u9st[:,:,itime].copy()
        H.v9[:,:] = v9st[:,:,itime].copy()
        H.w9[:,:] = w9st[:,:,itime].copy()
        H.t9[:,:] = t9st[:,:,itime].copy()
        H.p9[:,:] = p9st[:,:,itime].copy()
        H.qv9[:,:] = qv9st[:,:,itime].copy()
        H.ql9[:,:] = ql9st[:,:,itime].copy()
        H.u19[:,:] = u19st[:,:,itime].copy()
        H.v19[:,:] = v19st[:,:,itime].copy()
        H.w19[:,:] = w19st[:,:,itime].copy()
        H.t19[:,:] = t19st[:,:,itime].copy()
        H.p19[:,:] = p19st[:,:,itime].copy()
        H.qv19[:,:] = qv19st[:,:,itime].copy()
        H.ql19[:,:] = ql19st[:,:,itime].copy()

        hurr_adj()

    H.u[ 1:  ,1:-1]+=H.u1[()]  
    H.v[ 1:-1,1:-1]+=H.v1[()]  
    H.w[ 1:-1,1:  ]+=H.w1[()]  
    H.t[ 1:-1,1:-1]+=H.t1[()]  
    H.p[ 1:-1,1:-1]+=H.p1[()]  
    H.qv[1:-1,1:-1]+=H.qv1[()] 
    H.ql[1:-1,1:-1]+=H.ql1[()] 

    # limiting vertical levels
    H.u[:,top:]=0
    H.v[:,top:]=0
    H.w[:,top:]=0
    H.t[:,top:]=0
    H.p[:,top:]=0
    H.qv[:,top:]=0
    H.ql[:,top:]=0
    # end limiting

    dvarnd=[H.u , H.v , H.w , H.t , H.p , H.qv , H.ql]
    buf, dvar1d=nd2oned(dvarnd)

    return dvar1d

def grad_cost(avar1d):
    precon=get_precon()
    ax=avar1d*precon # convert from L-1*x to x

    obs_dicts=get_obs_dicts()
    res=0
    for obs_dict in obs_dicts: res+=grad_cost_obs(ax, obs_dict)
    res+=grad_cost_mod(ax)
    res*=precon # convert from dx to L*dx
    
    #print 'grad', np.linalg.norm(res)
    return res
