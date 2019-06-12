from peakondrake import *

dts = [0.01, 0.005, 0.001]
Ld = 40.
tmax = 500

for i, dt in enumerate(dts):
    code = 'skdv_dt_'+str(i)
    dt = dts[i]

    experiment(code, Ld, tmax,
               resolutions=200,
               dts=dt,
               sigmas=[0.002, 0.02, 0.2],
               seeds=range(20),
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
               field_ndump=int(tmax / (1 * dt)),
               allow_fail=True)
