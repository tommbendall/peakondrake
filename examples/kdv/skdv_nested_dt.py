from peakondrake import *

dts = [0.012, 0.01, 0.006, 0.004, 0.002]
code = 'skdv_nested_dt_'
Ld = 40.
tmax = 120


for i, dt in enumerate(dts):
    dt = dts[i]
    nXi_update = int(dt / np.min(dts))
    experiment(code+str(i), Ld, tmax,
               resolutions=200,
               dts=dt,
               sigmas=[0.002/sqrt(nXi_update), 0.02/sqrt(nXi_update), 0.2/sqrt(nXi_update)],
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
               field_ndump=int(tmax / (1 * dt)),
               allow_fail=True,
               nXi_updates=nXi_update)
