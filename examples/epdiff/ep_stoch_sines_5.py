from peakondrake import *


code = 'ep_peak_new_mu_sines_5'
Ld = 40.
dt = 0.0002
tmax = 40
i = 9
code += '_'+str(i)

experiment(code, Ld, tmax,
           resolutions=[200, 800, 1500, 4000, 10000],
           dts=dt,
           sigmas=0.5,
           seeds=range(10*i, 10*(i+1)),
           schemes='upwind',
           timesteppings='midpoint',
           ics='one_peak',
           num_Xis=6,
           Xi_family='sines',
           alphasq=1.0,
           c0=0.,
           gamma=0.,
           diagnostics=['mu'],
           fields_to_output=['uscalar', 'du'],
           ndump=int(tmax / (1000 * dt)),
           field_ndump=int(tmax / (100 * dt)),
           allow_fail=True)
