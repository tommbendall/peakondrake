from peakondrake.simulation import *
from firedrake import PeriodicIntervalMesh
from os import path, makedirs
from netCDF4 import Dataset
from collections import OrderedDict
import numpy as np
import types


def experiment(code, Ld, tmax, resolutions=[],
               dts=[], sigmas=[], seeds=[],
               schemes=[], timesteppings=[],
               alphasq=[], c0=[], gamma=[],
               ics=[], num_Xis=[], Xi_family=[],
               diagnostics=None, fields_to_output=None,
               ndump=-1, field_ndump=-1, allow_fail=False,
               t_kick=[], sigma_kick=0.0, smooth_t=None):

    # set up dumping
    dirname = 'results/'+code
    if path.exists(dirname):
        raise IOError("results directory '%s' already exists" % dirname)
    else:
        makedirs(dirname)

    expmt_dict = {'dt':(float, dts),
                  'resolution':(float, resolutions),
                  'seed':(int, seeds),
                  'sigma':(float, sigmas),
                  'scheme':('S1', schemes),
                  'timestepping':('S1', timesteppings),
                  'alphasq':(float, alphasq),
                  'c0':(float, c0),
                  'gamma':(float, gamma),
                  'ic':('S1', ics),
                  'num_Xis':(int, num_Xis),
                  'Xi_family':('S1', Xi_family)}

    # turns expmt_dict values into lists so they can be iterated over
    for key, value in expmt_dict.items():
        # need to not accidentally accept strings
        if not isinstance(value, str):
            try:
                iter(value[1])
            except TypeError:
                expmt_dict[key] = (value[0], [value[1]])
        else:
            expmt_dict[key] = (value[0], [value[1]])

    # we need a dictionary of parameters to form the loops over
    # also a dictionary of fixed parameters for all simulations
    # finally a dictionary of all parameters for individual simulations
    # the first should be an ordered dict
    variable_parameters = OrderedDict()
    fixed_parameters = {}
    simulation_parameters = {}

    # need resolutions to be first in variable parameters dict (to loop over first)
    if len(expmt_dict['resolution'][1]) > 1:
        variable_parameters['resolution'] = None

    # here we separate parameters into variable or fixed
    # this is done by asking if they are iterable, and being careful to ignore strings
    for key, value in expmt_dict.items():
        if not isinstance(value[1], str) and len(value[1]) > 1:
            variable_parameters[key] = value
            simulation_parameters[key] = None
        elif isinstance(value[1], str):
            fixed_parameters[key] = value
            simulation_parameters[key] = (value[1],)
        elif len(value[1]) == 1:
            fixed_parameters[key] = value
            simulation_parameters[key] = (value[1][0],)

    for key, value in zip(['dirname', 'Ld', 'tmax', 'ndump', 'field_ndump', 'allow_fail', 't_kick', 'sigma_kick', 'smooth_t'],
                          [dirname, Ld, tmax, ndump, field_ndump, allow_fail, t_kick, sigma_kick, smooth_t]):
        simulation_parameters[key] = (value,)

    output_arguments = ('time',) + tuple(variable_parameters.keys())

    for i, dt in enumerate(expmt_dict['dt'][1]):
        # make data file
        if len(expmt_dict['dt'][1]) > 1:
            file_name = dirname+'/data_dt'+str(i)+'.nc'
        else:
            file_name = dirname+'/data.nc'

        simulation_parameters['file_name'] = (file_name, )

        data_file = Dataset(file_name, 'w')

        # create dimensions and variables
        data_file.createDimension('time', None)
        data_file.createVariable('time', float, ('time',))
        data_file.createDimension('x', 100)
        data_file.createVariable('x', float, ('x',))
        data_file['x'][:] = np.linspace(0, Ld, num=100, endpoint=False)

        if 'resolution' in variable_parameters.keys():
            data_file.createDimension('resolution', len(variable_parameters['resolution'][1]))
            data_file.createVariable('deltax', float, ('resolution',))
            for j, res in enumerate(variable_parameters['resolution'][1]):
                data_file['deltax'][j:j+1] =  Ld / variable_parameters['resolution'][1][j]

        for key, values in variable_parameters.items():
            if key != 'resolution':
                data_file.createDimension(key, len(values[1]))
                data_file.createVariable(key, values[0], (key,))
                for j, value in enumerate(values[1]):
                    data_file[key][j:j+1] = value

        if diagnostics is not None:
            for output in diagnostics:
                if output == 'mu':
                    for i in range(10):
                        data_file.createVariable(output+'_'+str(i), float, output_arguments)
                elif output in ('a', 'b'):
                    data_file.createVariable(output, float, output_arguments+('x',))
                else:
                    data_file.createVariable(output, float, output_arguments)
        data_file.close()

        # we want to do a variable number of for loops so use recursive strategy

        N = len(variable_parameters)
        if N > 0:
            # measure length of loops
            # raise error if it's too long
            num_sims = np.prod([len(value[1]) for value in variable_parameters.values()])
            if num_sims > 1000:
                raise ValueError('You are asking to do %d simulations! Please reduce your loops.' % num_sims)

            # do the simulation loop
            if len(expmt_dict['resolution'][1]) == 1:
                resolution = simulation_parameters['resolution'][-1]
                mesh = PeriodicIntervalMesh(resolution, Ld)
                simulation_parameters['mesh'] = (mesh,)
            do_simulation_loop(N, variable_parameters, simulation_parameters,
                               diagnostics=diagnostics, fields_to_output=fields_to_output)
        else:
            resolution = simulation_parameters['resolution'][-1]
            mesh = PeriodicIntervalMesh(resolution, Ld)
            simulation_parameters['mesh'] = (mesh,)

            simulation(simulation_parameters, diagnostic_values=diagnostics,
                       fields_to_output=fields_to_output)

def do_simulation_loop(N, variable_parameters, simulation_parameters,
                       diagnostics=None, fields_to_output=None):
    """
    A recursive strategy for setting up a variable number of for loops
    for the experiment.

    :arg N: the level of loop.
    :arg variable_parameters: an OrderedDict of the parameters to be varied.
    :arg simulation_parameters: a dict storing the parameters for that simulation.
    :arg diagnostics: a list of diagnostic values to be output.
    :arg fields_to_output: a list of fields to be output.
    """

    # we do the loop in reverse order to get the resolution loop on the outside
    M = len(variable_parameters)
    key = list(variable_parameters.items())[M-N][0]

    # we must turn the ordered dict into a list to iterate through it
    for index, value in enumerate(list(variable_parameters.items())[M-N][1][1]):
        simulation_parameters[key] = (index, value)

        # make mesh if loop is a resolution loop
        if key == 'resolution':
            mesh = PeriodicIntervalMesh(value, simulation_parameters['Ld'][-1])
            simulation_parameters['mesh'] = (mesh,)

        # do recursion if we aren't finished yet
        if N > 1:
            do_simulation_loop(N-1, variable_parameters, simulation_parameters,
                               diagnostics=diagnostics, fields_to_output=fields_to_output)

        # finally do simulation
        elif N == 1:

            simulation(simulation_parameters,
                       diagnostic_values=diagnostics,
                       fields_to_output=fields_to_output)
