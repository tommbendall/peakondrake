from peakondrake import *

base_code = 'integrity_hydrodynamic_dt'
Ld = 40.
tmax = 50

dts = [0.00001, 0.0001, 0.001, 0.01]

for i, dt in enumerate(dts):

    code = base_code + '_'+str(i)

    experiment(code, Ld, tmax,
               resolutions=10000,
               dts=dt,
               sigmas=[0.0, 0.2],
               seeds=0,
               schemes='hydrodynamic',
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
               nXi_update=int(np.max(dts)/dt))
