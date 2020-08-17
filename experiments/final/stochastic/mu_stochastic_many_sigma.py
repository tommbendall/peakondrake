from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime
import sys

base_code = 'mu_stochastic_many_sigma'
Ld = 40.
tmax = 40
dt = 0.001
sigmas=[0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0]

if len(sys.argv) > 2:
    raise ValueError('We can only take one argument at the moment')

# sys.argv[0] is the name of the python file
i = sys.argv[1]
code = base_code+'_'+str(i)

starttime = datetime.now()

experiment(code, Ld, tmax,
           resolutions=2000,
           dts=dt,
           sigmas=sigmas[i],
           seeds=range(200),
           schemes='hydrodynamic',
           timesteppings='midpoint',
           ics='one_peak',
           num_Xis=1,
           Xi_family='constant',
           alphasq=1.0,
           c0=0.,
           gamma=0.,
           diagnostics=['mu', 'max_du', 'min_du'],
           fields_to_output=['du'],
           ndump=int(tmax / (2000 * dt)),
           field_ndump=int(tmax / (1 * dt)),
           allow_fail=True,
           nXi_updates=1)

endtime = datetime.now()
print(code)
print('Total runtime was %i' % (endtime - starttime).seconds)
