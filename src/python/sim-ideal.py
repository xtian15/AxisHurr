import numpy as np
import h5py as h5

import matplotlib.pyplot as plt
from subprocess import call
from datetime import datetime, timedelta

from axis_hurr import *

H=hurr_vars

hurr_in()
param_setup()
hurr_initial()

vmax=[]
etime=35000
for H.itime in range(1, etime+1):
    hurr_fwd()
    ws_max=np.sqrt(H.u9**2+H.v9**2).max()
    vmax.append(ws_max)
    
outname = 'Ideal-{:d}.h5'.format(etime)
with h5.File(outname, 'w') as outfile:
    outfile.create_dataset('u9', data=H.u9)
    outfile.create_dataset('v9', data=H.v9)
    outfile.create_dataset('w9', data=H.w9)
    outfile.create_dataset('t9', data=H.t9)
    outfile.create_dataset('p9', data=H.p9)
    outfile.create_dataset('qv9', data=H.qv9)
    outfile.create_dataset('ql9', data=H.ql9)
    outfile.create_dataset('vmax', data=vmax)
