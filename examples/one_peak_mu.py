from peakondrake import *


code = 'smooth_mu'
Ld = 40.
dt = 0.001
tmax = 50

experiment(code, Ld, tmax,
           resolutions=[500, 625, 750, 875, 1000],
           dts=dt,
           sigmas=0.0,
           seeds=0,
           schemes='upwind',
           timesteppings='midpoint',
           ics='one_peak',
           num_Xis=5,
           Xi_family='sines',
           alphasq=1.0,
           c0=0.,
           gamma=0.,
           diagnostics=['l2_m', 'max_jump_local', 'max_jump_global', 'max_du_loc', 'min_du_loc', 'max_du_smooth_loc', 'min_du_smooth_loc'],
           fields_to_output=['uscalar', 'du', 'jump_du', 'du_smooth'],
           ndump=int(tmax / (100 * dt)),
           field_ndump=int(tmax / (100 * dt)))
