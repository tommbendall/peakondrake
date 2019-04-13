from peakondrake import *

code = 'test_jump'
Ld = 40.
dt = 0.001
tmax = 0.01

experiment(code, Ld, tmax,
           resolutions=[100, 500, 1000, 5000, 10000],
           dts=dt,
           sigmas=0.0,
           seeds=0,
           schemes='upwind',
           timesteppings='midpoint',
           ics='one_peak',
           num_Xis=0,
           Xi_family='sines',
           alphasq=1.0,
           c0=0.,
           gamma=0.,
           diagnostics=['l2_m', 'max_jump_local', 'max_du_loc', 'min_du_loc'],
           fields_to_output=['uscalar', 'du'],
           ndump=int(tmax / (10 * dt)),
           field_ndump=int(tmax / (10 * dt)))
