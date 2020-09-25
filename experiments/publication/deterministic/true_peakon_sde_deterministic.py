from peakondrake import *
from netCDF4 import Dataset

base_code = 'true_deterministic_peakon_data'
Ld = 40.
tmax = 1.0
dt = 1e-6


experiment(base_code, Ld, tmax,
           resolutions=10000,
           dts=dt,
           sigmas=0.0,
           seeds=0,
           schemes='hydrodynamic',
           timesteppings='midpoint',
           ics='periodic_peakon',
           num_Xis=1,
           Xi_family='constant',
           alphasq=1.0,
           c0=0.,
           gamma=0.,
           diagnostics=['l2_u'],
           fields_to_output=[],
           ndump=int(tmax / (1 * dt)),
           field_ndump=int(tmax / (1 * dt)),
           allow_fail=True,
           peakon_equations=True,
           only_peakons=True,
           peakon_method='milstein',
           periodic=True)
