from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime

base_code = 'final_mu_strong_dt_sigma_20_convergence_'
Ld = 40.
tmax = 80
dts = [0.0005, 0.001]

starttime = datetime.now()

for i, dt in enumerate(dts):
    code = base_code+str(i)
    experiment(code, Ld, tmax,
               resolutions=1500,
               dts=dt,
               sigmas=0.2/sqrt(int(dt/np.min(dts))),
               seeds=0,
               schemes='hydrodynamic',
               timesteppings='midpoint',
               ics='one_peak',
               num_Xis=6,
               Xi_family='sines',
               alphasq=1.0,
               c0=0.,
               gamma=0.,
               diagnostics=['mu', 'l2_u'],
               fields_to_output=['du'],
               ndump=int(tmax / (2000 * dt)),
               field_ndump=int(tmax / (1 * dt)),
               allow_fail=True,
               nXi_updates=int(dt/np.min(dts)))

endtime = datetime.now()
print(code)
print('Total runtime was %i' % (endtime - starttime).seconds)
