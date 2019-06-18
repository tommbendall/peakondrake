from firedrake import FunctionSpace, Function, VectorFunctionSpace, Constant, SpatialCoordinate, as_vector, ConvergenceError
from collections import OrderedDict
from peakondrake.initial_conditions import *
from peakondrake.equations import *
from peakondrake.diagnostic_equations import *
from peakondrake.stochastic_functions import *
from peakondrake.outputting import *
from peakondrake.peakon_equations import *
import numpy as np

def simulation(simulation_parameters,
               diagnostic_values=None,
               fields_to_output=None):

    """
    This function sets up and runs a basic simulation of the stochastic
    Camassa-Holm equation.

    :arg simulation_parameters: a dictionary containing the parameters for
        the simulation. Each value should be a tuple, which will normally
        have only one element. Two-element tuples will also contain an index
        for outputting -- these are variables being looped over by the experiment.
    :arg diagnostic_values: a list of diagnostic numbers to output.
    :arg fields_to_output: a list of diagnostic fields to be dumped.
    """

    # ensure simulation parameters use ufl Constant
    for key in ['alphasq', 'c0', 'gamma']:
        if not isinstance(simulation_parameters[key][-1], Constant):
            if len(simulation_parameters[key]) == 2:
                simulation_parameters[key] = (simulation_parameters[key][0],
                                              Constant(simulation_parameters[key][1]))
            elif len(simulation_parameters[key]) == 1:
                simulation_parameters[key] = (Constant(simulation_parameters[key][-1]),)
            else:
                raise ValueError('The value of simulation_parameters %s seems to be of length %d. It should be either 1 or 2.' % (key, len(simulation_parameters[key])))

    mesh = simulation_parameters['mesh'][-1]
    dt = simulation_parameters['dt'][-1]
    tmax = simulation_parameters['tmax'][-1]
    ndump = simulation_parameters['ndump'][-1]
    field_ndump = simulation_parameters['field_ndump'][-1]
    file_name = simulation_parameters['file_name'][-1]
    scheme = simulation_parameters['scheme'][-1]
    allow_fail = simulation_parameters['allow_fail'][-1]
    only_peakons = simulation_parameters['only_peakons'][-1]

    allow_fail = True if only_peakons else allow_fail


    prognostic_variables = PrognosticVariables(scheme, mesh)
    diagnostic_variables = DiagnosticVariables(prognostic_variables, fields_to_output, diagnostic_values=diagnostic_values)

    build_initial_conditions(prognostic_variables, simulation_parameters)
    Xis = StochasticFunctions(prognostic_variables, simulation_parameters)

    equations = Equations(prognostic_variables, simulation_parameters)

    # set up peakon equations
    if simulation_parameters['peakon_equations'][-1] == True:
        peakon_equations = PeakonEquations(prognostic_variables, simulation_parameters)
    elif simulation_parameters['peakon_equations'][-1] == False:
        peakon_equations = None
    else:
        raise ValueError('peakon_equations must be True or False, not %s' % simulation_parameters['peakon_equations'][-1])

    outputting = Outputting(prognostic_variables,
                            diagnostic_variables,
                            simulation_parameters,
                            peakon_equations=peakon_equations,
                            diagnostic_values=diagnostic_values)

    diagnostic_equations = DiagnosticEquations(diagnostic_variables, prognostic_variables, outputting, simulation_parameters)
    diagnostic_equations.solve()



    dumpn = 0
    field_dumpn = 0
    # want the update to always happen on first time step (so results are consistent with earlier results)
    t = 0
    failed = True if only_peakons else False
    failed_time = np.nan

    # run simulation
    while (t < tmax - 0.5*dt):
        t += dt

        # update stochastic basis functions
        Xis.update(t)

        # solve problems
        if allow_fail:
            if not failed:
                try:
                    equations.solve()
                except ConvergenceError:
                    failed = True
                    failed_time = t
                    print("Solver failed at t = %f" % t)
        else:
            equations.solve()

        if peakon_equations is not None:
            peakon_equations.update()

        # output if necessary
        dumpn += 1
        field_dumpn += 1

        # do interpolations/projections for any outputing
        if dumpn == ndump or field_dumpn == field_ndump:
            if not failed:
                diagnostic_equations.solve()

        # output diagnostic values
        if dumpn == ndump:
            outputting.dump_diagnostics(t, failed=failed)
            dumpn -= ndump

        # output fields
        if field_dumpn == field_ndump:
            if not failed:
                outputting.dump_fields(t)
                field_dumpn -= field_ndump

    outputting.store_times(failed_time)

class PrognosticVariables(object):
    """
    Object setting up and storing the prognostic variables
    and their function spaces.

    :arg scheme: string specifying which equation set to use.
    :arg mesh: the mesh of the domain.
    """
    def __init__(self, scheme, mesh):

        self.scheme = scheme
        self.mesh = mesh
        self.fields = {}

        if scheme == 'upwind':
            self.Vm = FunctionSpace(mesh, "DG", 1)
            self.Vu = VectorFunctionSpace(mesh, "CG", 1)
            self.m = Function(self.Vm, name='m')
            self.u = Function(self.Vu, name='u')
            self.fields['u'] = self.u
            self.fields['m'] = self.m
        elif scheme == 'conforming':
            self.Vu = FunctionSpace(mesh, "CG", 1)
            self.u = Function(self.Vu, name='u')
            self.u0 = Function(self.Vu, name='u0')
            self.Vf = FunctionSpace(mesh, "CG", 1)
            self.fields['u'] = self.u
            self.Vm = FunctionSpace(mesh, "CG", 1)
            self.m = Function(self.Vm, name='m')
            self.fields['m'] = self.m
        elif scheme == 'hydrodynamic':
            self.Vu = FunctionSpace(mesh, "CG", 1)
            self.u = Function(self.Vu, name='u')
            self.Xi_x = Function(self.Vu, name='Xi_x')
            self.Xi_xx = Function(self.Vu, name='Xi_xx')
            self.fields['u'] = self.u
            self.fields['Xi_x'] = self.Xi_x
            self.fields['Xi_xx'] = self.Xi_xx
        elif scheme == 'test':
            self.Vu = FunctionSpace(mesh, "CG", 1)
            self.u = Function(self.Vu, name='u')
            self.fields['u'] = self.u
        else:
            raise ValueError('Scheme not recognised.')

        self.Xi = Function(self.Vu, name='Xi')
        self.fields['Xi'] = self.Xi
        self.pure_xi_list = []
        self.pure_xi_x_list = []

class DiagnosticVariables(object):
    """
    Object setting up and storing diagnostic variables.

    :arg prognostic_variables: a PrognosticVariables object.
    :arg fields_to_output: a list of diagnostic fields to use.
    """

    def __init__(self, prognostic_variables, fields_to_output, diagnostic_values=None):

        self.scheme = prognostic_variables.scheme
        self.mesh = prognostic_variables.mesh
        x, = SpatialCoordinate(self.mesh)
        self.coords = Function(FunctionSpace(self.mesh, "DG", 0)).interpolate(x)
        self.smooth_coords = Function(FunctionSpace(self.mesh, "CG", 5)).interpolate(x)
        self.fields = OrderedDict()
        self.dumpfields = OrderedDict()
        required_fields = []

        for field in fields_to_output:
            if field in ['uscalar', 'Xiscalar']:
                if self.scheme == 'upwind':
                    V = FunctionSpace(self.mesh, "CG", 1)
                else:
                    raise ValueError('Output field %s only usable with upwind scheme, not %s' % (field, self.scheme))
            elif field == 'du':
                V = FunctionSpace(self.mesh, "DG", 0)
            elif field == 'du_smooth':
                required_fields.extend(['du'])
                V = FunctionSpace(self.mesh, "CG", 5)
            elif field == 'jump_du':
                required_fields.extend(['du'])
                V = FunctionSpace(self.mesh, "CG", 1)
            elif field in ['F', 'u2_flux', 'a', 'b', 'kdv_1', 'kdv_2', 'kdv_3', 'm', 'u_xx', 'u_sde', 'u_sde_weak', 'u_sde_mean', 'u_sde_weak_mean']:
                V = FunctionSpace(self.mesh, "CG", 1)
            else:
                raise ValueError('Output field %s not recognised.' % field)

            self.fields[field] = Function(V, name=field)
            self.dumpfields[field] = self.fields[field]

        # if there is a diagnostic value that we want, we need to add any required fields
        # to the diagnostics so that they are calculated
        CG1 = FunctionSpace(self.mesh, "CG", 1)
        DG0 = FunctionSpace(self.mesh, "DG", 0)

        if diagnostic_values is not None:
            for diagnostic in diagnostic_values:
                if diagnostic == 'max_jump_local':
                    required_fields.extend(['jump_du', 'du'])
                elif diagnostic == 'max_jump_global':
                    required_fields.extend(['du'])
                elif diagnostic == 'max_du_loc':
                    required_fields.extend(['du'])
                elif diagnostic == 'min_du_loc':
                    required_fields.extend(['du'])
                elif diagnostic == 'max_du_smooth_loc':
                    required_fields.extend(['du', 'du_smooth'])
                elif diagnostic == 'min_du_smooth_loc':
                    required_fields.extend(['du', 'du_smooth'])
                elif diagnostic == 'mu':
                    required_fields.extend(['du'])
                elif diagnostic == 'a':
                    required_fields.extend(['a'])
                elif diagnostic == 'b':
                    required_fields.extend(['b'])
                elif diagnostic == 'l2_kdv_1':
                    required_fields.extend(['kdv_1'])
                elif diagnostic == 'l2_kdv_2':
                    required_fields.extend(['kdv_2'])
                elif diagnostic == 'l2_kdv_3':
                    required_fields.extend(['kdv_3'])
                elif diagnostic == 'mass_m':
                    if 'm' not in prognostic_variables.fields.keys():
                        required_fields.extend(['m'])
                elif diagnostic == 'm_max':
                    if 'm' not in prognostic_variables.fields.keys():
                        required_fields.extend(['m'])
                elif diagnostic == 'E_0':
                    required_fields.extend(['du'])
                elif diagnostic == 'E_2':
                    if 'm' not in prognostic_variables.fields.keys():
                        required_fields.extend(['m'])
                elif diagnostic == 'E_3':
                    required_fields.extend(['u_xx'])
                elif diagnostic == 'u_error_with_sde':
                    required_fields.extend(['u_sde'])
                elif diagnostic == 'u_error_weak':
                    required_fields.extend(['u_sde', 'u_sde_weak'])
                elif diagnostic == 'u_error_with_sde_mean':
                    required_fields.extend(['u_sde_mean'])
                elif diagnostic == 'u_error_weak_mean':
                    required_fields.extend(['u_sde_mean', 'u_sde_weak_mean'])

        for field in required_fields:
            if field not in self.fields.keys():
                if field == 'du':
                    self.fields[field] = Function(DG0)
                elif field == 'du_smooth':
                    self.fields[field] = Function(FunctionSpace(self.mesh, "CG", 5))
                else:
                    self.fields[field] = Function(CG1)
