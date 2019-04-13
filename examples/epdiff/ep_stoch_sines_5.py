from peakondrake import *


code = 'ep_peak_stochastic_sines_5'
Ld = 40.
dt = 0.001
tmax = 50

experiment(code, Ld, tmax,
           resolutions=[1000, 1250, 1500, 1750, 2000],
           dts=dt,
           sigmas=0.5,
           seeds=range(25),
           schemes='upwind',
           timesteppings='midpoint',
           ics='one_peak',
           num_Xis=6,
           Xi_family='sines',
           alphasq=1.0,
           c0=0.,
           gamma=0.,
           diagnostics=['l2_m', 'max_jump_local', 'max_jump_global', 'max_du_loc', 'min_du_loc'],
           fields_to_output=['uscalar', 'du', 'jump_du'],
           ndump=int(tmax / (100 * dt)),
           field_ndump=int(tmax / (100 * dt)))
