from firedrake import FunctionSpace, Function, VectorFunctionSpace, Constant
from peakondrake.initial_conditions import *
from peakondrake.equations import *
from peakondrake.diagnostic_equations import *
from peakondrake.stochastic_functions import *
from peakondrake.outputting import *

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


    prognostic_variables = PrognosticVariables(scheme, mesh)
    diagnostic_variables = DiagnosticVariables(prognostic_variables, fields_to_output, diagnostic_values=diagnostic_values)

    build_initial_conditions(prognostic_variables, simulation_parameters)
    Xis = StochasticFunctions(prognostic_variables, simulation_parameters)

    equations = Equations(prognostic_variables, simulation_parameters)
    diagnostic_equations = DiagnosticEquations(diagnostic_variables, prognostic_variables, simulation_parameters)
    diagnostic_equations.solve()

    outputting = Outputting(prognostic_variables,
                            diagnostic_variables,
                            simulation_parameters,
                            diagnostic_values=diagnostic_values)


    dumpn = 0
    field_dumpn = 0
    t = 0

    # run simulation
    while (t < tmax - 0.5*dt):
        t += dt

        # update stochastic basis functions
        Xis.update()

        # solve problems
        equations.solve()

        # output if necessary
        dumpn += 1
        field_dumpn += 1

        # do interpolations/projections for any outputing
        if dumpn == ndump or field_dumpn == field_ndump:
            diagnostic_equations.solve()

        # output diagnostic values
        if dumpn == ndump:
            outputting.dump_diagnostics(t)
            dumpn -= ndump

        # output fields
        if field_dumpn == field_ndump:
            outputting.dump_fields(t)
            field_dumpn -= field_ndump


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
            self.Vm = FunctionSpace(mesh, "CG", 1)
            self.Vu = FunctionSpace(mesh, "CG", 1)
            self.m = Function(self.Vm, name='m')
            self.u = Function(self.Vu, name='u')
            self.fields['u'] = self.u
            self.fields['m'] = self.m
        elif scheme == 'hydrodynamic':
            self.Vf = FunctionSpace(mesh, "CG", 1)
            self.Vu = FunctionSpace(mesh, "CG", 1)
            self.f = Function(self.Vf, name='f')
            self.u = Function(self.Vu, name='u')
            self.fields['u'] = self.u
            self.fields['f'] = self.f
        else:
            raise ValueError('Scheme not recognised.')

        self.Xi = Function(self.Vu, name='Xi')
        self.fields['Xi'] = self.Xi

class DiagnosticVariables(object):
    """
    Object setting up and storing diagnostic variables.

    :arg prognostic_variables: a PrognosticVariables object.
    :arg fields_to_output: a list of diagnostic fields to use.
    """

    def __init__(self, prognostic_variables, fields_to_output, diagnostic_values=None):

        self.scheme = prognostic_variables.scheme
        self.mesh = prognostic_variables.mesh
        self.fields = {}
        self.dumpfields = {}

        for field in fields_to_output:
            if field in ['uscalar', 'Xiscalar', 'jump_du']:
                if self.scheme == 'upwind':
                    V = FunctionSpace(self.mesh, "CG", 1)
                    self.fields[field] = Function(V, name=field)
                    self.dumpfields[field] = self.fields[field]
                else:
                    raise ValueError('Output field %s only usable with upwind scheme, not %s' % (field, self.scheme))
            elif field == 'du':
                V = FunctionSpace(self.mesh, "DG", 0)
                self.fields[field] = Function(V, name=field)
                self.dumpfields[field] = Function(V, name=field)

            else:
                raise ValueError('Output field %s not recognised.' % field)

        if diagnostic_values is not None:
            for diagnostic in diagnostic_values:
                if diagnostic == 'max_jump':
                    CG1 = FunctionSpace(self.mesh, "CG", 1)
                    DG0 = FunctionSpace(self.mesh, "DG", 0)
                    required_fields = ['jump_du', 'du']
                    for field in required_fields:
                        if field not in self.fields.keys():
                            if field == 'du':
                                self.fields[field] = Function(DG0)
                            else:
                                self.fields[field] = Function(CG1)
