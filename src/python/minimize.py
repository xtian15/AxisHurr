import numpy as np
import h5py as h5
import costfunction as cf
from subprocess import call
from matrix_reshape import *
from scipy.optimize import minimize, fmin_l_bfgs_b

fg=cf.get_fg()
shapes, avar1d=nd2oned(fg)
precon=cf.get_precon()
avar1d=avar1d/precon

res=fmin_l_bfgs_b(cf.cost, avar1d, fprime=cf.grad_cost, disp=True, 
                  pgtol=1.E-10)

#res=minimize(cf.cost, avar1d, method='CG', jac=cf.grad_cost, 
#             options={'disp': True})
print res
exit()

precon=cf.get_precon()
u9 , v9 , w9 , t9 , p9 , qv9 , ql9, \
u19, v19, w19, t19, p19, qv19, ql19=oned2nd(shapes, res[0]*precon)

with h5.File('outfile.h5', 'w') as outfile:
    u9=outfile.create_dataset('u9', data=u9)
    v9=outfile.create_dataset('v9', data=v9)
    w9=outfile.create_dataset('w9', data=w9)
    t9=outfile.create_dataset('t9', data=t9)
    p9=outfile.create_dataset('p9', data=p9)
    qv9=outfile.create_dataset('qv9', data=qv9)
    ql9=outfile.create_dataset('ql9', data=ql9)
