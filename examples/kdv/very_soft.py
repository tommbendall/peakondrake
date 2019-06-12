from peakondrake import *
from firedrake import sin, pi

dts = [0.012]
Ld = 40.
tmax = 1200

for i, dt in enumerate(dts):

    code = 'skdv_very_soft_more_output_dt_'+str(i)

    experiment(code, Ld, tmax,
               resolutions=200,
               dts=dt,
               sigmas=0.002,
               seeds=0,
               schemes='conforming',
               timesteppings='midpoint',
               ics='one_peak',
               num_Xis=3,
               Xi_family='sines',
               alphasq=0.0,
               c0=0.,
               gamma=1.0,
               diagnostics=['l2_u'],
               fields_to_output=[],
               ndump=int(tmax / (1000 * dt)),
               field_ndump=int(tmax / (1000 * dt)),
               allow_fail=True,
               nXi_update=(np.max(dts)/dt))
