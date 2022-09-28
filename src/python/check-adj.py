import numpy as np
import h5py as h5

from datetime import datetime, timedelta

# for debugging
#call(['make'])

from axis_hurr import *

H=hurr_vars

bfilename='Ideal-35000.h5'
with h5.File(bfilename, 'r') as bfile:
    bu9 = bfile['u9'][()]
    bv9 = bfile['v9'][()]
    bw9 = bfile['w9'][()]
    bt9 = bfile['t9'][()]
    bp9 = bfile['p9'][()]
    bqv9 = bfile['qv9'][()]
    bql9 = bfile['ql9'][()]

bu19=bu9[1:,1:-1]
bv19=bv9[1:-1,1:-1]
bw19=bw9[1:-1,1:]
bp19=bp9[1:-1,1:-1]
bt19=bt9[1:-1,1:-1]
bqv19=bqv9[1:-1,1:-1]
bql19=bql9[1:-1,1:-1]

n_tlmrun=300
u9st = np.zeros([bu9.shape[0], bu9.shape[1], n_tlmrun], float)
v9st = np.zeros([bv9.shape[0], bv9.shape[1], n_tlmrun], float)
w9st = np.zeros([bw9.shape[0], bw9.shape[1], n_tlmrun], float)
t9st = np.zeros([bt9.shape[0], bt9.shape[1], n_tlmrun], float)
p9st = np.zeros([bp9.shape[0], bp9.shape[1], n_tlmrun], float)
qv9st = np.zeros([bqv9.shape[0], bqv9.shape[1], n_tlmrun], float)
ql9st = np.zeros([bql9.shape[0], bql9.shape[1], n_tlmrun], float)
u19st = np.zeros([bu19.shape[0], bu19.shape[1], n_tlmrun], float)
v19st = np.zeros([bv19.shape[0], bv19.shape[1], n_tlmrun], float)
w19st = np.zeros([bw19.shape[0], bw19.shape[1], n_tlmrun], float)
t19st = np.zeros([bt19.shape[0], bt19.shape[1], n_tlmrun], float)
p19st = np.zeros([bp19.shape[0], bp19.shape[1], n_tlmrun], float)
qv19st = np.zeros([bqv19.shape[0], bqv19.shape[1], n_tlmrun], float)
ql19st = np.zeros([bql19.shape[0], bql19.shape[1], n_tlmrun], float)

for i in range(10):
    hurr_in()
    param_setup()
    hurr_initial()
    double_timer()
    
    pert=10**-i
    
    up = np.random.random_sample(H.u.shape)*bu9*pert
    vp = np.random.random_sample(H.v.shape)*bv9*pert
    wp = np.random.random_sample(H.w.shape)*bw9*pert
    tp = np.random.random_sample(H.t.shape)*bt9*pert
    pp = np.random.random_sample(H.p.shape)*bp9*pert
    qvp = np.random.random_sample(H.qv.shape)*bqv9*pert
    qlp = np.random.random_sample(H.ql.shape)*bql9*pert
    u1p = np.random.random_sample(H.u1.shape)*bu19*pert
    v1p = np.random.random_sample(H.v1.shape)*bv19*pert
    w1p = np.random.random_sample(H.w1.shape)*bw19*pert
    t1p = np.random.random_sample(H.t1.shape)*bt19*pert
    p1p = np.random.random_sample(H.p1.shape)*bp19*pert
    qv1p = np.random.random_sample(H.qv1.shape)*bqv19*pert
    ql1p = np.random.random_sample(H.ql1.shape)*bql19*pert

    H.u = up.copy()
    H.v = vp.copy()
    H.w = wp.copy()
    H.t = tp.copy()
    H.p = pp.copy()
    H.qv = qvp.copy()
    H.ql = qlp.copy()
    H.u1 = u1p.copy()
    H.v1 = v1p.copy()
    H.w1 = w1p.copy()
    H.t1 = t1p.copy()
    H.p1 = p1p.copy()
    H.qv1 = qv1p.copy()
    H.ql1 = ql1p.copy()

    H.u9 = bu9.copy()
    H.v9 = bv9.copy()
    H.w9 = bw9.copy()
    H.t9 = bt9.copy()
    H.p9 = bp9.copy()
    H.qv9 = bqv9.copy()
    H.ql9 = bql9.copy()
    
    H.u19 = bu19.copy()
    H.v19 = bv19.copy()
    H.w19 = bw19.copy()
    H.t19 = bt19.copy()
    H.p19 = bp19.copy()
    H.qv19 = bqv19.copy()
    H.ql19 = bql19.copy()

    for i in range(n_tlmrun):
        u9st[:,:,i]=H.u9[()].copy()
        v9st[:,:,i]=H.v9[()].copy()
        w9st[:,:,i]=H.w9[()].copy()
        t9st[:,:,i]=H.t9[()].copy()
        p9st[:,:,i]=H.p9[()].copy()
        qv9st[:,:,i]=H.qv9[()].copy()
        ql9st[:,:,i]=H.ql9[()].copy()
        u19st[:,:,i]=H.u19[()].copy()
        v19st[:,:,i]=H.v19[()].copy()
        w19st[:,:,i]=H.w19[()].copy()
        t19st[:,:,i]=H.t19[()].copy()
        p19st[:,:,i]=H.p19[()].copy()
        qv19st[:,:,i]=H.qv19[()].copy()
        ql19st[:,:,i]=H.ql19[()].copy()
        
        hurr_tlm()

    H.lhs=(H.u**2).sum()+\
        (H.v**2).sum()+\
        (H.w**2).sum()+\
        (H.t**2).sum()+\
        (H.p**2).sum()+\
        (H.qv**2).sum()+\
        (H.ql**2).sum()+\
        (H.u1**2).sum()+\
        (H.v1**2).sum()+\
        (H.w1**2).sum()+\
        (H.t1**2).sum()+\
        (H.p1**2).sum()+\
        (H.qv1**2).sum()+\
        (H.ql1**2).sum()
    
    for i in range(n_tlmrun-1,-1,-1):
        H.u9[()]=u9st[:,:,i]
        H.v9[()]=v9st[:,:,i]
        H.w9[()]=w9st[:,:,i]
        H.t9[()]=t9st[:,:,i]
        H.p9[()]=p9st[:,:,i]
        H.qv9[()]=qv9st[:,:,i]
        H.ql9[()]=ql9st[:,:,i]
        H.u19[()]=u19st[:,:,i]
        H.v19[()]=v19st[:,:,i]
        H.w19[()]=w19st[:,:,i]
        H.t19[()]=t19st[:,:,i]
        H.p19[()]=p19st[:,:,i]
        H.qv19[()]=qv19st[:,:,i]
        H.ql19[()]=ql19st[:,:,i]
        
        hurr_adj()
    

    H.rhs=(H.u*up).sum() +\
        (H.v*vp).sum() +\
        (H.w*wp).sum() +\
        (H.t*tp).sum() +\
        (H.p*pp).sum() +\
        (H.qv*qvp).sum() +\
        (H.ql*qlp).sum() +\
        (H.u1*u1p).sum() +\
        (H.v1*v1p).sum() +\
        (H.w1*w1p).sum() +\
        (H.t1*t1p).sum() +\
        (H.p1*p1p).sum() +\
        (H.qv1*qv1p).sum() +\
        (H.ql1*ql1p).sum()

    print(pert, H.lhs, H.rhs, (H.lhs-H.rhs)/H.lhs)
