import numpy as np
import h5py as h5

from datetime import datetime, timedelta
from axis_hurr import *
from Hoperator import *
from matrix_reshape import *
from costfunction import *
from obs_subs import *

bu9, bv9, bw9, bt9, bp9, bqv9, bql9 = get_fg()
# ----- analysis argument coming from an upstream subprogram -----
au9 = bu9.copy()
av9 = bv9.copy()
aw9 = bw9.copy()
at9 = bt9.copy()
ap9 = bp9.copy()
aqv9 = bqv9.copy()
aql9 = bql9.copy()
# ----- analysis argument coming from an upstream subprogram ----- end

up=np.zeros(au9.shape, float)
vp=np.zeros(av9.shape, float)
wp=np.zeros(aw9.shape, float)
tp=np.zeros(at9.shape, float)
pp=np.zeros(ap9.shape, float)
qvp=np.zeros(aqv9.shape, float)
qlp=np.zeros(aql9.shape, float)

avarnd=[au9 , av9 , aw9 , at9 , ap9 , aqv9 , aql9]

obs_dicts=get_obs_dicts()
obs_dict=obs_dicts[0]

shapes, avar1d=nd2oned(avarnd)
#Jo1=cost_obs(avar1d, obs_dict)
precon=get_precon()
avar1d/=precon
Jo1=cost(avar1d)

#dvar1d=grad_cost_obs(avar1d, obs_dict)
dvar1d=grad_cost(avar1d)
dvar1d/=precon
dvarnd=oned2nd(shapes, dvar1d)

res=[]
for i in range(1,17):
    np.random.seed(10)
    up[1:-1,bot:top]=bu9[1:-1,bot:top]*(10**-i)
    vp[1:-1,bot:top]=bv9[1:-1,bot:top]*(10**-i)
    wp[1:-1,bot:top]=bw9[1:-1,bot:top]*(10**-i)
    tp[1:-1,bot:top]=bt9[1:-1,bot:top]*(10**-i)
    pp[1:-1,bot:top]=bp9[1:-1,bot:top]*(10**-i)
    qvp[1:-1,bot:top]=bqv9[1:-1,bot:top]*(10**-i)
    qlp[1:-1,bot:top]=bql9[1:-1,bot:top]*(10**-i)

    perts=[up, vp, wp, tp, pp, qvp, qlp]

    pvarnd=[au9+up, av9+vp, aw9+wp, at9+tp, ap9+pp, aqv9+qvp, aql9+qlp]
    shapes, pvar1d=nd2oned(pvarnd)
    #Jo2=cost_obs(pvar1d, obs_dict)
    pvar1d/=precon
    Jo2=cost(pvar1d)

    deno=0
    for ivar in range(len(perts)):
        deno+=(dvarnd[ivar]*perts[ivar]).sum()
    nume=Jo2-Jo1

    print 10**-i, nume/deno
    res.append(nume/deno)
#print np.array(res)
exit()
with h5.File('res-gradcheck.h5', 'w') as outfile:
    dset=outfile.create_dataset('res', data=res)
