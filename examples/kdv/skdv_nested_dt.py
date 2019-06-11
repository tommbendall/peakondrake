from peakondrake import *

dts = [0.012, 0.006, 0.004, 0.003, 0.002]
i = 0
code = 'skdv_nested_dt_'+str(i)
Ld = 40.
dt = dts[i]
tmax = 120
nXi_update = int(0.012/dt)

experiment(code, Ld, tmax,
           resolutions=200,
           dts=dt,
           sigmas=0.2,
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
           allow_fail=True,
           nXi_update=nXi_update)
