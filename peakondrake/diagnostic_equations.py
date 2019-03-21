from firedrake import (Projector, Interpolator, as_vector, Constant,
                       VectorFunctionSpace, Function, dot)


class DiagnosticEquations(object):
    """
    An object setting up diagnostic equations and solvers for
    the stochastic Camassa-Holm equation.

    :arg diagnostic_variables: a DiagnosticVariables object.
    :arg prognostic_variables: a PrognosticVariables object.
    :arg simulation_parameters: a dictionary storing the simulation parameters.
    """

    def __init__(self, diagnostic_variables, prognostic_variables, simulation_parameters):

        self.diagnostic_variables = diagnostic_variables
        self.prognostic_variables = prognostic_variables
        self.simulation_parameters = simulation_parameters
        u = self.prognostic_variables.u
        ones = Function(VectorFunctionSpace(self.prognostic_variables.mesh, "CG", 1)).project(as_vector([Constant(1.0)]))
        self.interpolators = []
        self.projectors = []
        self.solvers = []

        for key, value in self.diagnostic_variables.fields.items():

            if key == 'uscalar':
                uscalar = self.diagnostic_variables.fields['uscalar']
                u_interpolator = Interpolator(dot(ones, u), uscalar)
                self.interpolators.append(u_interpolator)

            elif key == 'Xiscalar':
                Xi = self.prognostic_variables.Xi
                Xiscalar = self.diagnostic_variables.fields['Xiscalar']
                Xi_interpolator = Interpolator(dot(ones, Xi), Xiscalar)
                self.interpolators.append(Xi_interpolator)

            else:
                raise NotImplementedError('Diagnostic %s not yet implemented' % key)


    def solve(self):
        """
        Do interpolations, projections and solves.
        """
        for interpolator in self.interpolators:
            interpolator.interpolate()
        for projector in self.projectors:
            projector.project()
        for solver in self.solvers:
            solver.solve()
