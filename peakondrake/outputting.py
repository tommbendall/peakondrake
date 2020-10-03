from firedrake import (File, assemble, dx, dS, norm, dot,
                       Function, op2, errornorm, exp, sqrt)
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
        self.list_of_peakon_diagnostics = ['peakon_loc', 'peakon_min_du', 'peakon_max_du',
                                           'peakon_min_du_loc', 'peakon_max_du_loc',
                                           'peakon_max_u', 'peakon_mu', 'peakon_nu']


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

        self.out_string = ''
        for key, value in simulation_parameters.items():
            if len(value) > 1:
                self.out_string += str(key) + ' = %s, ' % str(value[0])

        # write initial wallclock time
        self.data_file['wallclock_time'][[slice(0,1)]+self.index_slices] = datetime.now().timestamp()

        # set up coordinate field
        if self.simulation_parameters['store_coordinates'][-1]:
            self.data_file['x'][:] = self.diagnostic_variables.coords.dat.data[:]


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
        elif 'm' in self.diagnostic_variables.fields.keys():
            m = self.diagnostic_variables.fields['m']

        alphasq = self.simulation_parameters['alphasq'][-1]
        Ld = self.simulation_parameters['Ld'][-1]
        periodic = self.simulation_parameters['periodic'][-1]
        # Need to evaluate periodic factor because we need a number from Firedrake constants
        periodic_factor = (1 - exp(-Ld/sqrt(alphasq)))
        periodic_factor = periodic_factor.evaluate(0,0,0,0)

        for diagnostic in self.diagnostic_values:
            if failed:
                output = [[np.nan], [np.nan]] if diagnostic == 'mu' else np.nan
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
            elif diagnostic == 'mass_u':
                output = assemble(u * dx)
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
            elif diagnostic == 'max_du':
                output = np.max(self.diagnostic_variables.fields['du'].dat.data[:])
            elif diagnostic == 'min_du':
                output = np.min(self.diagnostic_variables.fields['du'].dat.data[:])
            elif diagnostic == 'mu':
                old_mu, alt_mu = find_mus(u, self.diagnostic_variables.fields['du'], self.diagnostic_variables.coords)
                output = [old_mu, alt_mu]
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
                if periodic:
                    output = np.max(u.dat.data[:]) * periodic_factor
                else:
                    output = np.max(u.dat.data[:])
            elif diagnostic == 'q_pde':
                output = self.diagnostic_variables.coords.dat.data[np.argmax(u.dat.data[:])]
            elif diagnostic == 'm_max':
                output = np.max(m.dat.data[:])
            elif diagnostic == 'E_0':
                # use the calculated du
                du = self.diagnostic_variables.fields['du']
                output = assemble(0.5*(u**2 + alphasq*du**2)*dx)
            elif diagnostic == 'E_1':
                # just straightforwardly use u
                output = assemble(0.5*(u**2 + alphasq*u.dx(0)**2)*dx)
            elif diagnostic == 'E_2':
                # use m
                output = assemble(0.5*u*m*dx)
            elif diagnostic == 'E_3':
                # solve for uxx
                u_xx = self.diagnostic_variables.fields['u_xx']
                output = assemble(0.5*u*(u + alphasq*u_xx)*dx)
            elif diagnostic == 'u_error_with_sde':
                u_hat = self.diagnostic_variables.fields['u_sde']
                output = errornorm(u, u_hat)
            elif diagnostic == 'u_error_weak':
                u_hat_weak = self.diagnostic_variables.fields['u_sde_weak']
                output = norm(u_hat_weak)
            elif diagnostic == 'u_error_with_sde_mean':
                u_hat = self.diagnostic_variables.fields['u_sde_mean']
                output = errornorm(u, u_hat)
            elif diagnostic == 'u_error_weak_mean':
                u_hat_weak = self.diagnostic_variables.fields['u_sde_weak_mean']
                output = norm(u_hat_weak)
            elif diagnostic == 'u_field':
                output = u.dat.data[:]
            elif diagnostic == 'peakon_suite':
                output = peakon_diagnostics(u, self.diagnostic_variables.fields['du'], self.diagnostic_variables.coords)
            else:
                raise ValueError('Diagnostic %s not recgonised.' % diagnostic)

            if diagnostic in ('a', 'b'):
                self.data_file[diagnostic][[slice(self.t_idx,self.t_idx+1)]+self.index_slices] = output
            elif diagnostic == 'mu':
                # we cannot store arrays of mus, so have to do them each separately
                for i in range(4):
                    if i < len(output[0]):
                        # output[0] is the old mu, output[1] is the new mu
                        self.data_file[diagnostic+'_'+str(i)][[slice(self.t_idx,self.t_idx+1)]+self.index_slices] = output[0][i]
                        self.data_file['alt_'+diagnostic+'_'+str(i)][[slice(self.t_idx,self.t_idx+1)]+self.index_slices] = output[1][i]
                    else:
                        self.data_file[diagnostic+'_'+str(i)][[slice(self.t_idx,self.t_idx+1)]+self.index_slices] = np.nan
                        self.data_file['alt_'+diagnostic+'_'+str(i)][[slice(self.t_idx,self.t_idx+1)]+self.index_slices] = np.nan
            elif diagnostic == 'peakon_suite':
                for peakon_diag in self.list_of_peakon_diagnostics:
                    self.data_file[peakon_diag][[slice(self.t_idx,self.t_idx+1)]+self.index_slices] = output[peakon_diag]
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
    """
    In the alternative method, we return distance between peak and min df.
    In the original method, we return the distance between max df and min df.
    """

    list_of_mus = []
    list_of_xmin = []
    list_of_xmax = []
    alt_list_of_mus = []
    approx_peak_locations = []
    peak_locations = []
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
                peak_locations.append(0.5*(x.dat.data[i-1]+x.dat.data[i]))
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
            # only do this if one of the 4 highest peaks
            if len(approx_peak_locations) < 4 or peak_heights[i] > sorted_peak_heights[3]:
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


                if abs(xmin - xmax) < L / 2:
                    list_of_mus.append(abs(xmin - xmax))
                else:
                    list_of_mus.append(L - abs(xmin - xmax))
                if abs(xmin - peak_locations[i]) < L / 2:
                    alt_list_of_mus.append(abs(xmin - peak_locations[i]))
                else:
                    alt_list_of_mus.append(L - abs(xmin - peak_locations[i]))
                # these were just for sanity checks
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


        if len(peak_locations) > 0:
            if abs(xmin - peak_locations[0]) < L / 2:
                alt_list_of_mus.append(abs(xmin - peak_locations[0]))
            else:
                alt_list_of_mus.append(L - abs(xmin - peak_locations[0]))
        else:
            # if there is no recorded peak, take half of the old mu
            print('No peakon detected, using half-mu instead')
            if abs(xmin - xmax) < L / 2:
                alt_list_of_mus.append(abs(xmin - xmax)/2)
            else:
                alt_list_of_mus.append((L - abs(xmin - xmax))/2)

        list_of_xmax.append(xmax)
        list_of_xmin.append(xmin)

    return list_of_mus, alt_list_of_mus


def peakon_diagnostics(f, df, x):
    """
    This returns a set of diagnostics for determining whether the highest
    peak is a peakon and whether it is experiencing wave breaking.

    :arg x: the coordinates in DG0
    """

    L = np.max(x.dat.data[:]) + np.min(x.dat.data[:]) # x is DG0 coords so this gives the length

    dof_count = len(df.dat.data[:])
    search_length = int(dof_count/2)

    # find the location of the highest peak
    max_u = np.max(f.dat.data[:])
    peak_idx_CG1 = np.argmax(f.dat.data[:])  # index for CG1
    peak_idx_DG0_m = peak_idx_CG1 - 1
    peak_idx_DG0_p = peak_idx_CG1  # conveniently this should not exceed len(dat.data[:])

    # x is the coords in DG0
    # Need correct treatment if the peak is at the edge of the domain
    if x.dat.data[peak_idx_DG0_m] > x.dat.data[peak_idx_DG0_p]:
         # In this situation it should just be 0
        peak_loc = 0.5*(x.dat.data[peak_idx_DG0_m] + x.dat.data[peak_idx_DG0_p] - L)
    else:
        # Otherwise do normal way of finding a CG1 location from
        peak_loc = 0.5*(x.dat.data[peak_idx_DG0_m] + x.dat.data[peak_idx_DG0_p])

    # Find locations of local mins and maxs in du around the peak
    # Start with min du to left of peak and max du to right of peak
    guess_min_du = df.dat.data[peak_idx_DG0_m]
    guess_max_du = df.dat.data[peak_idx_DG0_p]

    min_du = np.nan
    max_du = np.nan
    x_min_du = np.nan
    x_max_du = np.nan

    # Find min du
    for idx in range(peak_idx_DG0_m+1, peak_idx_DG0_m+1+search_length):
        DG0_idx = idx - dof_count if (idx > dof_count - 1) else idx

        if df.dat.data[DG0_idx] > guess_min_du:
            min_du = guess_min_du
            min_du_idx = DG0_idx - 1
            x_min_du = x.dat.data[min_du_idx]

            break
        else:
            guess_min_du = df.dat.data[DG0_idx]

    # Find max du
    for idx in range(peak_idx_DG0_p-1, peak_idx_DG0_p-1-search_length, -1):
        DG0_idx = idx + dof_count if (idx < 0) else idx

        if df.dat.data[DG0_idx] < guess_max_du:
            max_du = guess_max_du
            max_du_idx = DG0_idx - 1
            x_max_du = x.dat.data[max_du_idx]

            break
        else:
            guess_max_du = df.dat.data[DG0_idx]

    # make final important diagnostics
    mu = (x_min_du - x_max_du) if (x_max_du < x_min_du) else (L - x_min_du + x_max_du)
    nu = (max_du - min_du) / mu

    diagnostic_dict = {'peakon_loc': peak_loc,
                       'peakon_min_du': min_du,
                       'peakon_max_du': max_du,
                       'peakon_min_du_loc': x_min_du,
                       'peakon_max_du_loc': x_max_du,
                       'peakon_max_u': max_u,
                       'peakon_mu': mu,
                       'peakon_nu': nu}

    return diagnostic_dict
