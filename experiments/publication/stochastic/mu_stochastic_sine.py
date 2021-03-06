from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime
import sys

base_code = 'new_mu_stochastic_sine'
Ld = 40.
tmax = 20
dt = 0.001

# sys.argv[0] is the name of the python file
i = int(sys.argv[1])
j = int(sys.argv[2])
index = 10*j + i
code = base_code+'_'+str(index)

starttime = datetime.now()

experiment(code, Ld, tmax,
           resolutions=[1000, 1500, 2000, 3000, 5000],
           dts=dt,
           sigmas=[0.05, 0.1, 0.2, 0.5, 1.0],
           seeds=range(5*index,5*(index+1)),
           schemes='hydrodynamic',
           timesteppings='midpoint',
           ics='proper_peak',
           num_Xis=6,
           Xi_family='sines',
           alphasq=1.0,
           c0=0.,
           gamma=0.,
           diagnostics=['mu', 'max_du', 'min_du'],
           fields_to_output=['du'],
           ndump=int(tmax / (1000 * dt)),
           field_ndump=int(tmax / (1 * dt)),
           allow_fail=True,
           nXi_updates=1)

endtime = datetime.now()
print(code)
print('Total runtime was %i' % (endtime - starttime).seconds)
