from peakondrake import *

code = 'test_conforming_ch'
Ld = 40.
dt = 0.001
tmax = 100

experiment(code, Ld, tmax,
           resolutions=200,
           dts=dt,
           sigmas=0.2,
           seeds=0,
           schemes='conforming',
           timesteppings='midpoint',
           ics='one_peak',
           num_Xis=1,
           Xi_family='sines',
           alphasq=1.0,
           c0=0.,
           gamma=0.0,
           diagnostics=['mu'],
           fields_to_output=['du'],
           ndump=int(tmax / (100 * dt)),
           field_ndump=int(tmax / (100 * dt)),
           allow_fail=True)
