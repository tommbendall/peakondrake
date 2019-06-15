from peakondrake import *


code = 'skdv_very_long'
Ld = 40.
dt = 0.001
tmax = 20000

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
           field_ndump=int(tmax / (1 * dt)),
           allow_fail=True)