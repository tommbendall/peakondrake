from peakondrake.simulation import *
from firedrake import PeriodicIntervalMesh, File
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
               ndump=-1, field_ndump=-1, nXi_updates=1, allow_fail=False,
               t_kick=[], sigma_kick=0.0, smooth_t=None,
               peakon_equations=False, only_peakons=False,
               true_peakon_data=None, true_mean_peakon_data=None,
               peakon_speed=None, expected_u=False, periodic=False,
               peak_width=1.0, peakon_method='ito_euler',
               fixed_dW=None):

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
                  'Xi_family':('S1', Xi_family),}

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

    # Stuff where there can be only one value that is fixed for all experiments
    for key, value in zip(['dirname', 'Ld', 'tmax', 'ndump', 'field_ndump', 'nXi_updates', 'allow_fail', 'smooth_t', 'peakon_equations', 'only_peakons', 'true_peakon_data', 'true_mean_peakon_data', 'peakon_speed', 'periodic', 'fixed_dW', 'peak_width', 'peakon_method'],
                          [dirname, Ld, tmax, ndump, field_ndump, nXi_updates, allow_fail, smooth_t, peakon_equations, only_peakons, true_peakon_data, true_mean_peakon_data, peakon_speed, periodic, fixed_dW, peak_width, peakon_method]):
        simulation_parameters[key] = (value,)

    output_arguments = ('time',) + tuple(variable_parameters.keys())
    list_of_peakon_diagnostics = ['peakon_loc', 'peakon_min_du', 'peakon_max_du',
                                       'peakon_min_du_loc', 'peakon_max_du_loc',
                                       'peakon_max_u', 'peakon_mu', 'peakon_nu']

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


        if 'resolution' in variable_parameters.keys():
            data_file.createDimension('resolution', len(variable_parameters['resolution'][1]))
            data_file.createVariable('deltax', float, ('resolution',))
            for j, res in enumerate(variable_parameters['resolution'][1]):
                data_file['deltax'][j:j+1] =  Ld / variable_parameters['resolution'][1][j]

            # create a fixed x dimension for a and b fields
            data_file.createDimension('x', 100)
            data_file.createVariable('x', float, ('x',))
            data_file['x'][:] = np.linspace(0, Ld, num=100, endpoint=False)
            simulation_parameters['store_coordinates'] = [False]

        else:
            # if the resolution is not varying, we can use the real resolution
            num_points = simulation_parameters['resolution'][-1]
            data_file.createDimension('x', num_points)
            data_file.createVariable('x', float, ('x',))
            simulation_parameters['store_coordinates'] = [True]

        for key, values in variable_parameters.items():
            if key != 'resolution':
                data_file.createDimension(key, len(values[1]))
                data_file.createVariable(key, values[0], (key,))
                for j, value in enumerate(values[1]):
                    data_file[key][j:j+1] = value

        # create diagnostic variables, with dimensions dependent upon parameters
        if diagnostics is not None:
            for output in diagnostics:
                if output == 'mu':
                    for i in range(4):
                        data_file.createVariable(output+'_'+str(i), float, output_arguments)
                        data_file.createVariable('alt_'+output+'_'+str(i), float, output_arguments)
                elif output in ('a', 'b', 'u_field'):
                    data_file.createVariable(output, float, output_arguments+('x',))
                elif output == 'peakon_suite':
                    for peakon_diag in list_of_peakon_diagnostics:
                        data_file.createVariable(peakon_diag, float, output_arguments)
                else:
                    data_file.createVariable(output, float, output_arguments)

        if peakon_equations:
            data_file.createVariable('p', float, output_arguments)
            data_file.createVariable('q', float, output_arguments)

        # create diagnostics for wallclock time and failure time
        data_file.createDimension('wallclock_times', 2)
        data_file.createVariable('wallclock_time', float, ('wallclock_times',) + tuple(variable_parameters.keys()))
        data_file.createVariable('failed_time', float, tuple(variable_parameters.keys()))

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
                               diagnostics=diagnostics, fields_to_output=fields_to_output,
                               expected_u=expected_u)
        else:
            resolution = simulation_parameters['resolution'][-1]
            mesh = PeriodicIntervalMesh(resolution, Ld)
            simulation_parameters['mesh'] = (mesh,)

            simulation(simulation_parameters, diagnostic_values=diagnostics,
                       fields_to_output=fields_to_output)

def do_simulation_loop(N, variable_parameters, simulation_parameters,
                       diagnostics=None, fields_to_output=None, expected_u=False):
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
    have_setup = False

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
                               diagnostics=diagnostics, fields_to_output=fields_to_output, expected_u=expected_u)

        # finally do simulation
        elif N == 1:

            if expected_u:
                this_u = simulation(simulation_parameters,
                           diagnostic_values=diagnostics,
                           fields_to_output=fields_to_output,
                           expected_u=True)
                if have_setup:
                    Eu.assign(counter * Eu + this_u)
                    counter.assign(counter + 1)
                    Eu.assign(Eu / counter)
                else:
                    scheme = simulation_parameters['scheme'][-1]
                    mesh = simulation_parameters['mesh'][-1]
                    prognostic_variables = PrognosticVariables(scheme, mesh)
                    Eu = Function(prognostic_variables.Vu, name='expected u').assign(this_u)
                    counter = Constant(1.0)
                    have_setup = True
            else:
                simulation(simulation_parameters,
                           diagnostic_values=diagnostics,
                           fields_to_output=fields_to_output)


    if expected_u:
        expected_u_file = File(simulation_parameters['dirname'][-1]+'/expected_u.pvd')
        expected_u_file.write(Eu)
