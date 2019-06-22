from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime
from os import path, makedirs

code = 'peakon_speed'
Ld = 40.
tmax = 1000
dt = 0.01
dxs = [400, 800, 1000, 1500, 2000, 5000, 10000]

def unfurl_q(q_data, Ld):
    new_q = [q_data[0]]
    laps = 0
    for q in q_data[1:]:
        if new_q[-1] > 9*Ld/10 + laps * Ld and q < Ld / 10:
            laps += 1
        new_q.append(q + laps * Ld)

    return np.array(new_q)

starttime = datetime.now()

experiment(code, Ld, tmax,
           resolutions=dxs,
           dts=dt,
           sigmas=0.0,
           seeds=0,
           schemes='hydrodynamic',
           timesteppings='midpoint',
           ics='peakon',
           num_Xis=0,
           Xi_family='sines',
           alphasq=1.0,
           c0=0.,
           gamma=0.,
           diagnostics=['p_pde', 'q_pde'],
           fields_to_output=[],
           ndump=int(tmax / (1000 * dt)),
           field_ndump=int(tmax / (1 * dt)),
           allow_fail=True,
           peakon_equations=False)

endtime = datetime.now()
print(code)
print('Total runtime was %i' % (endtime - starttime).seconds)


data = Dataset('results/'+code+'/data.nc', 'r')
speed_data = Dataset('results/'+code+'/speed_data.nc', 'w')
speed_data.createDimension('dx', len(dxs))
speed_data.createVariable('dx', float, ('dx',))
speed_data.createVariable('speed', float, ('dx',))

final_positions = []
speeds = []
for i, dx in enumerate(dxs):
    q = data['q_pde'][:,i]
    unfurled_q = unfurl_q(q, Ld)
    speed_data['dx'][i:i+1] = dx
    speed_data['speed'][i:i+1] = (unfurled_q[-1] + Ld / 2) / tmax
