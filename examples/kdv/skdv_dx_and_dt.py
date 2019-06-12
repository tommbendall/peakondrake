from peakondrake import *

dxs = [200, 800, 2000]
dts = [0.01, 0.0025, 0.001]
Ld = 40.
tmax = 500


for i, (dt, dx) in enumerate(zip(dts, dxs)):
    code = 'skdv_dx_and_dt_'+str(i)

    experiment(code, Ld, tmax,
               resolutions=dx,
               dts=dt,
               sigmas=[0.002, 0.02, 0.2],
               seeds=range(20),
               schemes='conforming',
               timesteppings='midpoint',
               ics='one_peak',
               num_Xis=1,
               Xi_family='gaussians',
               alphasq=0.0,
               c0=0.,
               gamma=1.0,
               diagnostics=['l2_u'],
               fields_to_output=[],
               ndump=int(tmax / (1000 * dt)),
               field_ndump=int(tmax / (1 * dt)),
               allow_fail=True)
