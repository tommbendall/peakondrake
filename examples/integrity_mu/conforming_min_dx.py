from peakondrake import *

code = 'integrity_conforming_min_dx'
Ld = 40.
tmax = 50
dt = 0.00001

experiment(code, Ld, tmax,
           resolutions=[1000, 2000, 5000, 10000],
           dts=dt,
           sigmas=[0.0, 0.2],
           seeds=0,
           schemes='conforming',
           timesteppings='midpoint',
           ics='one_peak',
           num_Xis=3,
           Xi_family='sines',
           alphasq=1.0,
           c0=0.,
           gamma=0.,
           diagnostics=['mu'],
           fields_to_output=['du'],
           ndump=int(tmax / (1000 * dt)),
           field_ndump=int(tmax / (10 * dt)),
           allow_fail=True,
           nXi_update=1)
