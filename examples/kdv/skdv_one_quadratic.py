from peakondrake import *


code = 'skdv_one_quadratic'
Ld = 40.
dt = 0.001
tmax = 500

experiment(code, Ld, tmax,
           resolutions=200,
           dts=dt,
           sigmas=0.02,
           seeds=range(50),
           schemes='conforming',
           timesteppings='midpoint',
           ics='one_peak',
           num_Xis=1,
           Xi_family='quadratic',
           alphasq=0.0,
           c0=0.,
           gamma=1.0,
           diagnostics=['l2_u'],
           fields_to_output=[],
           ndump=int(tmax / (100 * dt)),
           field_ndump=int(tmax / (100 * dt)),
           allow_fail=True)
