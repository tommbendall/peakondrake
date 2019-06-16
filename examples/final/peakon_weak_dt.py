from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime

base_code = 'final_peakon_convergence_weak_dt'
Ld = 40.
tmax = 200
dts = [0.01, 0.001, 0.0001]

starttime = datetime.now()

for i, dt in enumerate(dts):
    code = base_code+'_'+str(i)
    experiment(code, Ld, tmax,
               resolutions=800,
               dts=dt,
               sigmas=0.02,
               seeds=range(20),
               schemes='hydrodynamic',
               timesteppings='midpoint',
               ics='peakon',
               num_Xis=1,
               Xi_family='sines',
               alphasq=1.0,
               c0=0.,
               gamma=0.,
               diagnostics=['p_pde', 'q_pde'],
               fields_to_output=[],
               ndump=int(tmax / (2000 * dt)),
               field_ndump=int(tmax / (1 * dt)),
               allow_fail=True,
               peakon_equations=True)

endtime = datetime.now()
print(code)
print('Total runtime was %i' % (endtime - starttime).seconds)
