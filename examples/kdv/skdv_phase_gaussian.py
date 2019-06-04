from peakondrake import *


code = 'skdv_phase_gaussian'
Ld = 40.
dt = 0.001
tmax = 500

experiment(code, Ld, tmax,
           resolutions=200,
           dts=dt,
           sigmas=0.2,
           seeds=range(50),
           schemes='conforming',
           timesteppings='midpoint',
           ics='one_peak',
           num_Xis=1,
           Xi_family='gaussians',
           alphasq=0.0,
           c0=0.,
           gamma=1.0,
           diagnostics=['l2_u', 'a', 'b'],
           fields_to_output=[],
           ndump=int(tmax / (10000 * dt)),
           field_ndump=int(tmax / (10000 * dt)),
           allow_fail=True)