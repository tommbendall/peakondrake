from peakondrake import *

code = 'test_one_run'
Ld = 40.
dt = 0.001
tmax = 0.01

experiment(code, Ld, tmax,
           resolutions=20,
           dts=dt,
           sigmas=0.25,
           seeds=0,
           schemes='upwind',
           timesteppings='ssprk3',
           ics='one_peak',
           num_Xis=7,
           Xi_family='sines',
           alphasq=1.0,
           c0=0.,
           gamma=0.,
           diagnostics=['h1_u', 'h1_m', 'l2_u', 'l2_m', 'energy', 'max_jump'],
           fields_to_output=['uscalar', 'Xiscalar'],
           ndump=int(tmax / (10 * dt)),
           field_ndump=int(tmax / (10 * dt)))
