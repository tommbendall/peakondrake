from peakondrake import *

dts = [0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001, 0.00005]
i = 0
code = 'skdv_dt_'+str(i)
Ld = 40.
dt = dts[i]
tmax = 500

experiment(code, Ld, tmax,
           resolutions=200,
           dts=dt,
           sigmas=0.2,
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
