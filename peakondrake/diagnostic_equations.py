from firedrake import (Projector, Interpolator, as_vector, Constant,
                       VectorFunctionSpace, Function, dot, TestFunction,
                       TrialFunction, LinearVariationalProblem, dx,
                       LinearVariationalSolver, dS, jump, VectorElement,
                       NonlinearVariationalSolver, NonlinearVariationalProblem,
                       FunctionSpace, SpatialCoordinate, sqrt, split, MixedFunctionSpace,
                       TestFunctions, conditional, exp)
from netCDF4 import Dataset


class DiagnosticEquations(object):
    """
    An object setting up diagnostic equations and solvers for
    the stochastic Camassa-Holm equation.

    :arg diagnostic_variables: a DiagnosticVariables object.
    :arg prognostic_variables: a PrognosticVariables object.
    :arg outputting: an Outputting object.
    :arg simulation_parameters: a dictionary storing the simulation parameters.
    """

    def __init__(self, diagnostic_variables, prognostic_variables, outputting, simulation_parameters):

        self.diagnostic_variables = diagnostic_variables
        self.prognostic_variables = prognostic_variables
        self.outputting = outputting
        self.simulation_parameters = simulation_parameters
        Dt = Constant(simulation_parameters['dt'][-1])
        Ld = simulation_parameters['Ld'][-1]
        u = self.prognostic_variables.u
        Xi = self.prognostic_variables.Xi
        Vu = u.function_space()
        ones = Function(VectorFunctionSpace(self.prognostic_variables.mesh, "CG", 1)).project(as_vector([Constant(1.0)]))
        self.to_update_constants = False
        self.interpolators = []
        self.projectors = []
        self.solvers = []

        # do peakon data checks here
        true_peakon_data = simulation_parameters['true_peakon_data'][-1]
        if true_peakon_data is not None:
            self.true_peakon_file = Dataset('results/'+true_peakon_data+'/data.nc', 'r')
            # check length of file is correct
            ndump = simulation_parameters['ndump'][-1]
            tmax = simulation_parameters['tmax'][-1]
            dt = simulation_parameters['dt'][-1]
            if len(self.true_peakon_file['time'][:]) != int(tmax/(ndump*dt)):
                raise ValueError('If reading in true peakon data, the dump frequency must be the same as that used for the true peakon data.')
            if self.true_peakon_file['p'][:].shape != (int(tmax/(ndump*dt)),):
                raise ValueError('True peakon data must have same shape as proposed data!')

        # do peakon data checks here
        true_mean_peakon_data = simulation_parameters['true_mean_peakon_data'][-1]
        if true_mean_peakon_data is not None:
            self.true_mean_peakon_file = Dataset('results/'+true_mean_peakon_data+'/data.nc', 'r')
            # check length of file is correct
            ndump = simulation_parameters['ndump'][-1]
            tmax = simulation_parameters['tmax'][-1]
            dt = simulation_parameters['dt'][-1]
            if len(self.true_mean_peakon_file['time'][:]) != int(tmax/(ndump*dt)):
                raise ValueError('If reading in true peakon data, the dump frequency must be the same as that used for the true peakon data.')
            if self.true_mean_peakon_file['p'][:].shape != (int(tmax/(ndump*dt)),):
                raise ValueError('True peakon data must have same shape as proposed data!')

        mesh = u.function_space().mesh()
        x, = SpatialCoordinate(mesh)
        alphasq = simulation_parameters['alphasq'][-1]

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
                gamma = simulation_parameters['gamma'][-1]
                u2_flux = self.diagnostic_variables.fields['u2_flux']
                xis = self.prognostic_variables.pure_xi_list
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

                flux_expr = 0.0*x
                for xi, xi_x, xi_xxx in zip(xis, xis_x, xis_xxx):
                    flux_expr += (6*u.dx(0)*xi + 12*u*xi_x + gamma*xi_xxx) * (6*u.dx(0)*xi + 24*u*xi_x + gamma*xi_xxx)
                projector = Projector(flux_expr, u2_flux)
                self.projectors.append(projector)

            elif key == 'a':
                # find  6 * u_x * Xi + gamma * Xi_xxx
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
                # find 12 * u * Xi_x
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

            elif key == 'kdv_1':
                # find the first part of the kdv form
                u0 = prognostic_variables.u0
                uh = (u + u0) / 2
                us = Dt * uh + sqrt(Dt) * Xi
                psi = TestFunction(Vu)
                du_1 = self.diagnostic_variables.fields['kdv_1']

                eqn = psi * du_1 * dx - 6 * psi.dx(0) * uh * us * dx
                prob = NonlinearVariationalProblem(eqn, du_1)
                solver = NonlinearVariationalSolver(prob)
                self.solvers.append(solver)

            elif key == 'kdv_2':
                # find the second part of the kdv form
                u0 = prognostic_variables.u0
                uh = (u + u0) / 2
                us = Dt * uh + sqrt(Dt) * Xi
                psi = TestFunction(Vu)
                du_2 = self.diagnostic_variables.fields['kdv_2']

                eqn = psi * du_2 * dx + 6 * psi * uh * us.dx(0) * dx
                prob = NonlinearVariationalProblem(eqn, du_2)
                solver = NonlinearVariationalSolver(prob)
                self.solvers.append(solver)

            elif key == 'kdv_3':
                # find the third part of the kdv form
                u0 = prognostic_variables.u0
                uh = (u + u0) / 2
                us = Dt * uh + sqrt(Dt) * Xi
                du_3 = self.diagnostic_variables.fields['kdv_3']
                gamma = simulation_parameters['gamma'][-1]

                phi = TestFunction(Vu)
                F = Function(Vu)

                eqn = (phi * F * dx + phi.dx(0) * us.dx(0) * dx)
                prob = NonlinearVariationalProblem(eqn, F)
                solver = NonlinearVariationalSolver(prob)
                self.solvers.append(solver)

                self.projectors.append(Projector(-gamma * F.dx(0), du_3))

                # nu = TestFunction(Vu)
                # back_eqn = nu * du_3 * dx - gamma * nu.dx(0) * F * dx
                # back_prob = NonlinearVariationalProblem(back_eqn, du_3)
                # back_solver = NonlinearVariationalSolver(back_prob)
                # self.solvers.append(solver)

            elif key == 'm':

                m = self.diagnostic_variables.fields['m']
                phi = TestFunction(Vu)
                eqn = phi * m * dx - phi * u * dx - alphasq * phi.dx(0) * u.dx(0) * dx
                prob = NonlinearVariationalProblem(eqn, m)
                solver = NonlinearVariationalSolver(prob)
                self.solvers.append(solver)

            elif key == 'u_xx':

                u_xx = self.diagnostic_variables.fields['u_xx']
                phi = TestFunction(Vu)
                eqn = phi * u_xx * dx + phi.dx(0) * u_xx.dx(0) * dx
                prob = NonlinearVariationalProblem(eqn, u_xx)
                solver = NonlinearVariationalSolver(prob)
                self.solvers.append(solver)

            elif key == 'u_sde':
                self.to_update_constants = True
                self.p = Constant(1.0)
                self.q = Constant(Ld/2)

                u_sde = self.diagnostic_variables.fields['u_sde']
                expr = conditional(x < self.q - Ld / 2, self.p * exp(-(x-self.q+Ld)/sqrt(alphasq)),
                                   conditional(x < self.q + Ld / 2, self.p * exp(-sqrt((self.q-x)**2/alphasq)),
                                               self.p * exp(-(self.q+Ld-x)/sqrt(alphasq))))
                self.interpolators.append(Interpolator(expr, u_sde))

            elif key == 'u_sde_weak':
                u_sde = self.diagnostic_variables.fields['u_sde']
                u_sde_weak = self.diagnostic_variables.fields['u_sde_weak']
                psi = TestFunction(Vu)

                eqn = psi * u_sde_weak * dx - psi * (u - u_sde) * dx
                prob = NonlinearVariationalProblem(eqn, u_sde_weak)
                solver = NonlinearVariationalSolver(prob)
                self.solvers.append(solver)

            elif key == 'u_sde_mean':
                self.to_update_constants = True
                self.p = Constant(1.0)
                self.q = Constant(Ld/2)

                u_sde = self.diagnostic_variables.fields['u_sde_mean']
                expr = conditional(x < self.q - Ld / 2, self.p * exp(-(x-self.q+Ld)/sqrt(alphasq)),
                                   conditional(x < self.q + Ld / 2, self.p * exp(-sqrt((self.q-x)**2/alphasq)),
                                               self.p * exp(-(self.q+Ld-x)/sqrt(alphasq))))
                self.interpolators.append(Interpolator(expr, u_sde))

            elif key == 'u_sde_weak_mean':
                u_sde = self.diagnostic_variables.fields['u_sde_mean']
                u_sde_weak = self.diagnostic_variables.fields['u_sde_weak_mean']
                psi = TestFunction(Vu)

                eqn = psi * u_sde_weak * dx - psi * (u - u_sde) * dx
                prob = NonlinearVariationalProblem(eqn, u_sde_weak)
                solver = NonlinearVariationalSolver(prob)
                self.solvers.append(solver)

            else:
                raise NotImplementedError('Diagnostic %s not yet implemented' % key)

    def update_constants(self):
        """
        Update p and q constants from true peakon data.
        """
        self.p.assign(self.true_peakon_file['p'][self.outputting.t_idx])
        self.q.assign(self.true_peakon_file['q'][self.outputting.t_idx])

    def solve(self):
        """
        Do interpolations, projections and solves.
        """

        if self.to_update_constants:
            self.update_constants()
        for interpolator in self.interpolators:
            interpolator.interpolate()
        for projector in self.projectors:
            projector.project()
        for solver in self.solvers:
            solver.solve()
