from peakondrake import *


code = 'ep_peak_stochastic_sines_10'
Ld = 40.
dt = 0.0002
tmax = 40

experiment(code, Ld, tmax,
           resolutions=[200, 800, 1500, 4000, 10000],
           dts=dt,
           sigmas=1.0,
           seeds=range(10),
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
           ndump=int(tmax / (1000 * dt)),
           field_ndump=int(tmax / (100 * dt)),
           allow_fail=True)
