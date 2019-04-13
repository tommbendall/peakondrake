from firedrake import (File, assemble, dx, dS, norm, dot,
                       Function, op2)
from netCDF4 import Dataset
import numpy as np

class Outputting(object):
    """
    An object for outputting diagnostics and fields.

    :arg prognostic_variables: a PrognosticVariables object.
    :arg diagnostic_variables: a DiagnosticVariables object.
    :arg simulation_parameters: a dictionary containing the simulation parameters.
    :arg diagnostic_values: a list of diagnostics to use.
    """

    def __init__(self, prognostic_variables,
                 diagnostic_variables,
                 simulation_parameters,
                 diagnostic_values=None):

        self.ndump = simulation_parameters['ndump'][-1]
        self.field_ndump = simulation_parameters['field_ndump'][-1]
        self.prognostic_variables = prognostic_variables
        self.diagnostic_variables = diagnostic_variables
        self.simulation_parameters = simulation_parameters
        file_name = simulation_parameters['file_name'][-1]
        dirname = simulation_parameters['dirname'][-1]
        self.data_file = Dataset(file_name, 'a')

        # set up things for dumping diagnostics
        if diagnostic_values is not None:
            if isinstance(diagnostic_values, str):
                dimensions = self.data_file[diagnostic_values].dimensions
            else:
                dimensions = self.data_file[diagnostic_values[0]].dimensions

        index_list = []
        variable_list = []
        for dimension in dimensions:
            if dimension != 'time':
                variable_list.append(dimension)
                index_list.append(simulation_parameters[dimension][0])

        self.index_slices = [slice(index, index+1) for index in index_list]
        self.t_idx = 0
        self.diagnostic_values = diagnostic_values


        # set up things for dumping fields
        if self.field_ndump > 0:
            field_file_name = dirname+'/fields'
            for dimension, index in zip(variable_list, index_list):
                field_file_name += '_'+str(dimension)+str(index)
            field_file_name += '_output.pvd'
            self.field_file = File(field_file_name)

            prognostic_variables = [value for value in self.prognostic_variables.fields.values()]
            diagnostic_variables = [value for value in self.diagnostic_variables.dumpfields.values()]
            self.dumpfields = prognostic_variables + diagnostic_variables

            self.field_file.write(*self.dumpfields, t=0)

        self.out_string = ''
        for key, value in simulation_parameters.items():
            if len(value) > 1:
                self.out_string += str(key) + ' = %s, ' % str(value[0])

    def dump_diagnostics(self, t):
        """
        Dump the diagnostic values.

        :arg t: time.
        """

        print(self.out_string, 't = %.3f' % t)
        self.data_file.variables['time'][self.t_idx:self.t_idx + 1] = t

        u = self.prognostic_variables.u
        if 'm' in self.prognostic_variables.fields.keys():
            m = self.prognostic_variables.m

        alphasq = self.simulation_parameters['alphasq'][-1]

        for diagnostic in self.diagnostic_values:
            if diagnostic == 'energy':
                output = assemble((dot(u, u) + alphasq*dot(u.dx(0),u.dx(0)))*dx)
            elif diagnostic == 'l2_m':
                output = norm(m, norm_type='L2')
            elif diagnostic == 'l2_u':
                output = norm(u, norm_type='L2')
            elif diagnostic == 'h1_u':
                output = norm(u, norm_type='H1')
            elif diagnostic == 'h1_m':
                output = norm(m, norm_type='H1')
            elif diagnostic == 'mass_m':
                output = assemble(m * dx)
            elif diagnostic == 'mass_m2':
                output = assemble(m * m * dx)
            elif diagnostic == 'min_u':
                output = 1 if norm(u_min, norm_type='L2') > 1e-10 else 0
            elif diagnostic == 'max_jump_local':
                output = find_max(self.diagnostic_variables.fields['jump_du'], self.diagnostic_variables.coords)[0]
            elif diagnostic == 'max_jump_global':
                output = find_max(self.diagnostic_variables.fields['du'], self.diagnostic_variables.coords)[0] - find_min(self.diagnostic_variables.fields['du'], self.diagnostic_variables.coords)[0]
            elif diagnostic == 'max_du_loc':
                output = find_max(self.diagnostic_variables.fields['du'], self.diagnostic_variables.coords)[1]
            elif diagnostic == 'min_du_loc':
                output = find_min(self.diagnostic_variables.fields['du'], self.diagnostic_variables.coords)[1]
            elif diagnostic == 'max_du_smooth_loc':
                output = find_max(self.diagnostic_variables.fields['du_smooth'], self.diagnostic_variables.smooth_coords)[1]
            elif diagnostic == 'min_du_smooth_loc':
                output = find_min(self.diagnostic_variables.fields['du_smooth'], self.diagnostic_variables.smooth_coords)[1]
            else:
                raise ValueError('Diagnostic %s not recgonised.' % diagnostic)

            self.data_file[diagnostic][[slice(self.t_idx,self.t_idx+1)]+self.index_slices] = output

        self.t_idx += 1



    def dump_fields(self, t):
        """
        Dump the diagnostic fields.

        :arg t: time.
        """

        if self.field_ndump > 0:
            self.field_file.write(*self.dumpfields, time=t)


def find_min(f, x):
    fmin = np.min(f.dat.data[:])
    min_idx = np.argmin(f.dat.data[:])
    xmin = x.dat.data[min_idx]
    return (fmin, xmin)

def find_max(f, x):
    fmax = np.max(f.dat.data[:])
    max_idx = np.argmax(f.dat.data[:])
    xmax = x.dat.data[max_idx]
    return (fmax, xmax)
