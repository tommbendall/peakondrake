from peakondrake import *


code = 'deterministic_one_peak_jump'
Ld = 40.
dt = 0.0001
tmax = 50

experiment(code, Ld, tmax,
           resolutions=[200, 400, 800, 1600, 3200, 6400, 12800],
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
           diagnostics=['l2_m', 'max_jump'],
           fields_to_output=['uscalar', 'du', 'jump_du'],
           ndump=int(tmax / (100 * dt)),
           field_ndump=int(tmax / (100 * dt)))
