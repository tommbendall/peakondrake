from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime
from os import path, makedirs

Ld = 40.
tmax = 10
dt_true = 1e-2
dts = [1e-1, 2e-1, 4e-1]
seeds = range(5)
nout = 25

base_true_code = 'true_weak'
base_pde_code = 'final_peakon_convergence_weak_dt'

starttime = datetime.now()

for seed in seeds:
    true_code = base_true_code+'_'+str(seed)

    experiment(true_code, Ld, tmax,
               resolutions=10000,
               dts=dt_true,
               sigmas=0.02,
               seeds=0,
               schemes='hydrodynamic',
               timesteppings='midpoint',
               ics='peakon',
               num_Xis=1,
               Xi_family='double_sines',
               alphasq=1.0,
               c0=0.,
               gamma=0.,
               diagnostics=['l2_u'],
               fields_to_output=[],
               ndump=int(tmax / (nout * dt_true)),
               field_ndump=int(tmax / (1 * dt_true)),
               allow_fail=True,
               peakon_equations=True,
               only_peakons=True)


nt = nout
true_mean_code = 'true_weak_mean'

# set up dumping
dirname = 'results/'+true_mean_code
if path.exists(dirname):
    raise IOError("results directory '%s' already exists" % dirname)
else:
    makedirs(dirname)

mean_data = Dataset('results/'+true_mean_code+'/data.nc', 'w')
mean_data.createDimension('time', nt)
mean_data.createVariable('time', float, ('time',))
mean_data.createVariable('p', float, ('time',))
mean_data.createVariable('q', float, ('time',))


mean_p = np.zeros(nout)
mean_q = np.zeros(nout)
for seed in seeds:
    true_code = base_true_code+'_'+str(seed)
    data = Dataset('results/'+true_code+'/data.nc', 'r')
    mean_p += data['p'][:]
    mean_q += data['q'][:]
mean_p[:] = mean_p[:]/len(seeds)
mean_q[:] = mean_q[:]/len(seeds)
mean_data['time'][:] = data['time'][:]


for seed in seeds:
    true_code = base_true_code+'_'+str(seed)
    for i, dt in enumerate(dts):
        code = base_pde_code+'_'+str(seed)+'_'+str(i)
        experiment(code, Ld, tmax,
                   resolutions=1000,
                   dts=dt,
                   sigmas=0.02/sqrt(int(dt/dt_true)),
                   seeds=seed,
                   schemes='hydrodynamic',
                   timesteppings='midpoint',
                   ics='peakon',
                   num_Xis=1,
                   Xi_family='sines',
                   alphasq=1.0,
                   c0=0.,
                   gamma=0.,
                   diagnostics=['p_pde', 'q_pde', 'u_error_with_sde', 'u_error_weak', 'u_error_with_sde_mean', 'u_error_weak_mean'],
                   fields_to_output=['u_sde', 'u_sde_weak', 'u_sde_mean', 'u_sde_weak_mean'],
                   ndump=int(tmax / (nout * dt)),
                   field_ndump=int(tmax / (4 * dt)),
                   allow_fail=True,
                   peakon_equations=False,
                   true_peakon_data=true_code,
                   true_mean_peakon_data=true_mean_code,
                   nXi_updates=int(dt/dt_true))

endtime = datetime.now()
print(code)
print('Total runtime was %i' % (endtime - starttime).seconds)
