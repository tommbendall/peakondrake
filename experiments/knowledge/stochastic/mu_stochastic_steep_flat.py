from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime
import sys

base_code = 'mu_stochastic_steep_flat'
Ld = 40.
tmax = 20
dt = 0.0005

# sys.argv[0] is the name of the python file
i = int(sys.argv[1])
j = int(sys.argv[2])
index = 10*j + i
code = base_code+'_'+str(index)

starttime = datetime.now()

experiment(code, Ld, tmax,
           resolutions=[1000, 1500, 2000, 2500, 3000],
           dts=dt,
           sigmas=[0.05, 0.2, 0.5],
           seeds=range(5*index,5*(index+1)),
           schemes='hydrodynamic',
           timesteppings='midpoint',
           ics='new_peak',
           num_Xis=1,
           Xi_family='constant',
           alphasq=1.0,
           c0=0.,
           gamma=0.,
           diagnostics=['max_du', 'min_du', 'peakon_suite'],
           fields_to_output=['du'],
           ndump=int(tmax / (1000 * dt)),
           field_ndump=int(tmax / (1 * dt)),
           allow_fail=True,
           nXi_updates=1,
           peak_width=(1/6))

endtime = datetime.now()
print(code)
print('Total runtime was %i' % (endtime - starttime).seconds)
