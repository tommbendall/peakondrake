from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime

base_code = 'final_peakon_convergence_strong_dt'
Ld = 40.
tmax = 200
dts = [0.0001, 0.0002, 0.0004]

starttime = datetime.now()

for i, dt in enumerate(dts):
    code = base_code+'_'+str(i)
    experiment(code, Ld, tmax,
               resolutions=750,
               dts=dt,
               sigmas=0.02,
               seeds=range(0),
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
               peakon_equations=True,
               nXi_updates=int(dt/np.min(dts)))

endtime = datetime.now()
print(code)
print('Total runtime was %i' % (endtime - starttime).seconds)
