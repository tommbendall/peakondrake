from firedrake import (File, assemble, dx, dS, norm, dot,
                       Function, op2)
from netCDF4 import Dataset
import numpy as np
from datetime import datetime

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
                 peakon_equations=None,
                 diagnostic_values=None):

        self.ndump = simulation_parameters['ndump'][-1]
        self.field_ndump = simulation_parameters['field_ndump'][-1]
        self.prognostic_variables = prognostic_variables
        self.diagnostic_variables = diagnostic_variables
        self.simulation_parameters = simulation_parameters
        self.peakon_equations = peakon_equations
        file_name = simulation_parameters['file_name'][-1]
        dirname = simulation_parameters['dirname'][-1]
        self.data_file = Dataset(file_name, 'a')

        # set up things for dumping diagnostics
        if diagnostic_values is not None:
            if isinstance(diagnostic_values, str):
                dimensions = self.data_file[diagnostic_values].dimensions
            else:
                if diagnostic_values[0] != 'mu':
                    dimensions = self.data_file[diagnostic_values[0]].dimensions
                else:
                    dimensions = self.data_file[diagnostic_values[0]+'_'+str(0)].dimensions

        index_list = []
        variable_list = []
        for dimension in dimensions:
            if dimension not in ('time', 'x'):
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

        # write initial wallclock time
        self.data_file['wallclock_time'][[slice(0,1)]+self.index_slices] = datetime.now().timestamp()

    def dump_diagnostics(self, t, failed=False):
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
            if failed:
                output = [np.nan] if diagnostic == 'mu' else np.nan
            elif diagnostic == 'energy':
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
            elif diagnostic == 'mu':
                output = find_mus(u, self.diagnostic_variables.fields['du'], self.diagnostic_variables.coords)
            elif diagnostic == 'a':
                output = self.diagnostic_variables.fields['a'].at(self.data_file['x'][:], tolerance=1e-6)
            elif diagnostic == 'b':
                output = self.diagnostic_variables.fields['b'].at(self.data_file['x'][:], tolerance=1e-6)
            elif diagnostic == 'l2_kdv_1':
                output = norm(self.diagnostic_variables.fields['kdv_1'], norm_type='L2')
            elif diagnostic == 'l2_kdv_2':
                output = norm(self.diagnostic_variables.fields['kdv_2'], norm_type='L2')
            elif diagnostic == 'l2_kdv_3':
                output = norm(self.diagnostic_variables.fields['kdv_3'], norm_type='L2')
            elif diagnostic == 'p_pde':
                output = np.max(u.dat.data[:])
            elif diagnostic == 'q_pde':
                output = self.diagnostic_variables.coords.dat.data[np.argmax(u.dat.data[:])]
            else:
                raise ValueError('Diagnostic %s not recgonised.' % diagnostic)

            if diagnostic in ('a', 'b'):
                self.data_file[diagnostic][[slice(self.t_idx,self.t_idx+1)]+self.index_slices] = output
            elif diagnostic == 'mu':
                # we cannot store arrays of mus, so have to do them each separately
                for i, mu in enumerate(output):
                    self.data_file[diagnostic+'_'+str(i)][[slice(self.t_idx,self.t_idx+1)]+self.index_slices] = mu
            else:
                self.data_file[diagnostic][[slice(self.t_idx,self.t_idx+1)]+self.index_slices] = output

        if self.peakon_equations is not None:
            self.data_file['p'][[slice(self.t_idx,self.t_idx+1)]+self.index_slices] = self.peakon_equations.p
            self.data_file['q'][[slice(self.t_idx,self.t_idx+1)]+self.index_slices] = self.peakon_equations.q

        self.t_idx += 1



    def dump_fields(self, t):
        """
        Dump the diagnostic fields.

        :arg t: time.
        """

        if self.field_ndump > 0:
            self.field_file.write(*self.dumpfields, time=t)

    def store_times(self, failed_time):
        """
        Store the final wallclock time and the time at which the run failed.
        """

        self.data_file['wallclock_time'][[slice(1,2)]+self.index_slices] = datetime.now().timestamp()
        if len(self.index_slices) == 0:
            self.data_file['failed_time'][0] = failed_time
        else:
            self.data_file['failed_time'][self.index_slices] = failed_time




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

def find_mus(f, df, x):
    list_of_mus = []
    list_of_xmin = []
    list_of_xmax = []
    approx_peak_locations = []
    approx_peak_indices = []
    peak_heights = []
    L = np.max(x.dat.data[:])

    dof_count = len(df.dat.data[:])
    # count number of peaks
    for i in range(dof_count):
        if df.dat.data[i-1] > 0 and df.dat.data[i] < 0:
            # check that this is actually a peak
            if f.at([0.5*(x.dat.data[i-1]+x.dat.data[i])], tolerance=1e-06) > 2*np.mean(f.dat.data[:]):
                approx_peak_locations.append(x.dat.data[i])
                approx_peak_indices.append(i)
                peak_heights.append(f.at([0.5*(x.dat.data[i-1]+x.dat.data[i])], tolerance=1e-06))

    sorted_peak_heights = peak_heights.copy()
    sorted_peak_heights.sort(reverse=True)

    # divide domain into different peak zones
    dividing_points = []
    dividing_indices = []
    if len(approx_peak_locations) > 1:
        for i in range(len(approx_peak_locations)):
            if approx_peak_locations[i] > approx_peak_locations[i-1]:
                dividing_indices.append(int(np.floor(0.5*(approx_peak_indices[i-1]+approx_peak_indices[i]))))
            else:
                dividing_indices.append(int(np.floor(0.5*((approx_peak_indices[i-1]+approx_peak_indices[i]+dof_count))))
                                        % dof_count)
            dividing_points.append(x.dat.data[dividing_indices[-1]])

    # now look for max and min of df, in slices defined by dividing points
    if len(approx_peak_locations) > 1:
        for i in range(len(approx_peak_locations)):
            if dividing_indices[i] > dividing_indices[i-1]:
                max_idx = np.argmax(df.dat.data[dividing_indices[i-1]:dividing_indices[i]]) + dividing_indices[i-1]
                min_idx = np.argmin(df.dat.data[dividing_indices[i-1]:dividing_indices[i]]) + dividing_indices[i-1]
            else:
                # if dividing_indices[i] == 0 we can't simply slice
                if dividing_indices[i] == 0:
                    max_idx = np.argmax(df.dat.data[dividing_indices[i-1]:]) + dividing_indices[i-1]
                    min_idx = np.argmin(df.dat.data[dividing_indices[i-1]:]) + dividing_indices[i-1]
                # similarly with dividing_indices[i-1] == len(approx_peak_locations) - 1
                elif dividing_indices[i-1] == len(approx_peak_locations) - 1:
                    max_idx = np.argmax(df.dat.data[0:dividing_indices[i-1]]) + dividing_indices[i-1]
                    min_idx = np.argmin(df.dat.data[0:dividing_indices[i-1]]) + dividing_indices[i-1]
                # otherwise, our zone is at the start and end of the list, so we need to compare each patch of the list
                else:
                    max_idx_lower = np.argmax(df.dat.data[0:dividing_indices[i]])
                    max_idx_upper = np.argmax(df.dat.data[dividing_indices[i-1]:]) + dividing_indices[i-1]
                    min_idx_lower = np.argmin(df.dat.data[0:dividing_indices[i]])
                    min_idx_upper = np.argmin(df.dat.data[dividing_indices[i-1]:]) + dividing_indices[i-1]
                    max_idx = max_idx_lower if df.dat.data[max_idx_lower] > df.dat.data[max_idx_upper] else max_idx_upper
                    min_idx = min_idx_lower if df.dat.data[min_idx_lower] < df.dat.data[min_idx_upper] else min_idx_upper
            xmax = x.dat.data[max_idx]
            xmin = x.dat.data[min_idx]

            # only append if one of the 10 highest peaks
            if len(approx_peak_locations) < 11 or peak_heights[i] > sorted_peak_heights[10]:
                if abs(xmin - xmax) < L / 2:
                    list_of_mus.append(abs(xmin - xmax))
                else:
                    list_of_mus.append(L - abs(xmin - xmax))
                list_of_xmax.append(xmax)
                list_of_xmin.append(xmin)
    else:
        # just do straightforward thing if there is only one peak
        max_idx = np.argmax(df.dat.data[:])
        xmax = x.dat.data[max_idx]
        min_idx = np.argmin(df.dat.data[:])
        xmin = x.dat.data[min_idx]
        if abs(xmin - xmax) < L / 2:
            list_of_mus.append(abs(xmin - xmax))
        else:
            list_of_mus.append(L - abs(xmin - xmax))
        list_of_xmax.append(xmax)
        list_of_xmin.append(xmin)

    # print(list_of_mus, list_of_xmax, approx_peak_locations, list_of_xmin)

    return list_of_mus
