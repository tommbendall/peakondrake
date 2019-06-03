from firedrake import (Projector, Interpolator, as_vector, Constant,
                       VectorFunctionSpace, Function, dot, TestFunction,
                       TrialFunction, LinearVariationalProblem, dx,
                       LinearVariationalSolver, dS, jump, VectorElement,
                       NonlinearVariationalSolver, NonlinearVariationalProblem,
                       FunctionSpace, SpatialCoordinate)


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

            elif key == 'du':
                if type(u.function_space().ufl_element()) == VectorElement:
                    u_to_project = self.diagnostic_variables.fields['uscalar']
                else:
                    u_to_project = u
                du = self.diagnostic_variables.fields['du']
                du_projector = Projector(u_to_project.dx(0), du)
                self.projectors.append(du_projector)

            elif key == 'jump_du':
                du = self.diagnostic_variables.fields['du']
                jump_du = self.diagnostic_variables.fields['jump_du']
                V = jump_du.function_space()
                jtrial = TrialFunction(V)
                psi = TestFunction(V)
                Lj = psi('+') * abs(jump(du)) * dS
                aj = psi('+') * jtrial('+') * dS
                jprob = LinearVariationalProblem(aj, Lj, jump_du)
                jsolver = LinearVariationalSolver(jprob)
                self.solvers.append(jsolver)

            elif key == 'du_smooth':
                du = self.diagnostic_variables.fields['du']
                du_smooth = self.diagnostic_variables.fields['du_smooth']
                projector = Projector(du, du_smooth)
                self.projectors.append(projector)

            elif key == 'u2_flux':
                mesh = u.function_space().mesh()
                gamma = simulation_parameters['gamma'][-1]
                u2_flux = self.diagnostic_variables.fields['u2_flux']
                xis = self.prognostic_variables.pure_xis
                xis_x = []
                xis_xxx = []
                CG1 = FunctionSpace(mesh, "CG", 1)
                psi = TestFunction(CG1)
                for xi in xis:
                    xis_x.append(Function(CG1).project(xi.dx(0)))
                for xi_x in xis_x:
                    xi_xxx = Function(CG1)
                    form = (psi * xi_xxx + psi.dx(0) * xi_x.dx(0)) * dx
                    prob = NonlinearVariationalProblem(form, xi_xxx)
                    solver = NonlinearVariationalSolver(prob)
                    solver.solve()
                    xis_xxx.append(xi_xxx)

                x, = SpatialCoordinate(mesh)
                flux_expr = 0.0*x
                for xi, xi_x, xi_xxx in zip(xis, xis_x, xis_xxx):
                    flux_expr += (6*u.dx(0)*xi + 12*u*xi_x + gamma*xi_xxx) * (6*u.dx(0)*xi + 24*u*xi_x + gamma*xi_xxx)
                projector = Projector(flux_expr, u2_flux)
                self.projectors.append(projector)

            elif key == 'a':
                mesh = u.function_space().mesh()
                gamma = simulation_parameters['gamma'][-1]
                a_flux = self.diagnostic_variables.fields['a']
                xis = self.prognostic_variables.pure_xis
                xis_x = []
                xis_xxx = []
                CG1 = FunctionSpace(mesh, "CG", 1)
                psi = TestFunction(CG1)
                for xi in xis:
                    xis_x.append(Function(CG1).project(xi.dx(0)))
                for xi_x in xis_x:
                    xi_xxx = Function(CG1)
                    form = (psi * xi_xxx + psi.dx(0) * xi_x.dx(0)) * dx
                    prob = NonlinearVariationalProblem(form, xi_xxx)
                    solver = NonlinearVariationalSolver(prob)
                    solver.solve()
                    xis_xxx.append(xi_xxx)

                x, = SpatialCoordinate(mesh)
                a_expr = 0.0*x
                for xi, xi_x, xi_xxx in zip(xis, xis_x, xis_xxx):
                    a_expr += 6*u.dx(0)*xi + gamma*xi_xxx
                projector = Projector(a_expr, a_flux)
                self.projectors.append(projector)

            elif key == 'b':
                mesh = u.function_space().mesh()
                gamma = simulation_parameters['gamma'][-1]
                b_flux = self.diagnostic_variables.fields['b']
                xis = self.prognostic_variables.pure_xis

                x, = SpatialCoordinate(mesh)
                b_expr = 0.0*x
                for xi, xi_x, xi_xxx in zip(xis, xis_x, xis_xxx):
                    b_expr += 12*u*xi.dx(0)
                projector = Projector(b_expr, b_flux)
                self.projectors.append(projector)

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
