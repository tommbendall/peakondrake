from firedrake import (Interpolator, Constant, as_vector, sin,
                       cos, exp, FunctionSpace, VectorFunctionSpace,
                       pi, SpatialCoordinate, Function, conditional,
                       dx, TestFunction, NonlinearVariationalSolver,
                       NonlinearVariationalProblem)
import numpy as np

class StochasticFunctions(object):
    """
    The stochastic functions object.

    :arg prognostic_variables: a PrognosticVariables object.
    :arg simulation_parameters: a dictionary of the simulation parameters.
    """

    def __init__(self, prognostic_variables, simulation_parameters):

        mesh = simulation_parameters['mesh'][-1]
        x, = SpatialCoordinate(mesh)
        Ld = simulation_parameters['Ld'][-1]
        self.scheme = simulation_parameters['scheme'][-1]

        self.dt = simulation_parameters['dt'][-1]
        self.num_Xis = simulation_parameters['num_Xis'][-1]
        self.Xi_family = simulation_parameters['Xi_family'][-1]
        self.dXi = prognostic_variables.dXi
        self.dWs = [Constant(0.0) for dw in range(self.num_Xis)]
        self.dW_nums = prognostic_variables.dW_nums
        self.Xi_functions = []
        self.nXi_updates = simulation_parameters['nXi_updates'][-1]
        self.smooth_t = simulation_parameters['smooth_t'][-1]
        self.fixed_dW = simulation_parameters['fixed_dW'][-1]

        if self.smooth_t is not None and self.nXi_updates > 1:
            raise ValueError('Prescribing forcing and including multiple Xi updates are not compatible.')

        if self.smooth_t is not None or self.fixed_dW is not None:
            print('WARNING: Remember to change sigma to sigma * sqrt(dt) with the prescribed forcing option or the fixed_dW option.')


        seed = simulation_parameters['seed'][-1]
        np.random.seed(seed)

        # make sure sigma is a Constant
        if self.num_Xis != 0:
            if isinstance(simulation_parameters['sigma'][-1], Constant):
                self.sigma = simulation_parameters['sigma'][-1]
            else:
                self.sigma = Constant(simulation_parameters['sigma'][-1])
        else:
            self.sigma = Constant(0.0)

        self.pure_xi_list = prognostic_variables.pure_xi_list
        self.pure_xi_x_list = prognostic_variables.pure_xi_x_list
        self.pure_xi_xx_list = prognostic_variables.pure_xi_xx_list
        self.pure_xi_xxx_list = prognostic_variables.pure_xi_xxx_list
        self.pure_xi_xxxx_list = prognostic_variables.pure_xi_xxxx_list
        for xi in range(self.num_Xis):
            self.pure_xi_list.append(Function(self.dXi.function_space()))
            self.pure_xi_x_list.append(Function(self.dXi.function_space()))
            self.pure_xi_xx_list.append(Function(self.dXi.function_space()))
            self.pure_xi_xxx_list.append(Function(self.dXi.function_space()))
            self.pure_xi_xxxx_list.append(Function(self.dXi.function_space()))


        if self.Xi_family == 'sines':
            for n in range(self.num_Xis):
                if (n+1) % 2 == 1:
                    self.Xi_functions.append(self.sigma * sin(2*(n+1)*pi*x/Ld))
                else:
                    self.Xi_functions.append(self.sigma * cos(2*(n+1)*pi*x/Ld))

        elif self.Xi_family == 'double_sines':
            for n in range(self.num_Xis):
                if (n+1) % 2 == 1:
                    self.Xi_functions.append(self.sigma * sin(4*(n+1)*pi*x/Ld))
                else:
                    self.Xi_functions.append(self.sigma * cos(4*(n+1)*pi*x/Ld))

        elif self.Xi_family == 'high_freq_sines':
            for n in range(self.num_Xis):
                if (n+1) % 2 == 1:
                    self.Xi_functions.append(self.sigma * sin((2*(n+1)+10)*pi*x/Ld))
                else:
                    self.Xi_functions.append(self.sigma * cos((2*(n+1)+10)*pi*x/Ld))

        elif self.Xi_family == 'gaussians':
            for n in range(self.num_Xis):
                self.Xi_functions.append(self.sigma * 0.5*self.num_Xis*exp(-((x-Ld*(n+1)/(self.num_Xis +1.0))/2.)**2))

        elif self.Xi_family == 'quadratic':
            if self.num_Xis > 1:
                raise NotImplementedError('Quadratic Xi not yet implemented for more than one Xi')
            else:
                self.Xi_functions.append(32/(Ld*Ld)*conditional(x > Ld/4,
                                                     conditional(x > 3*Ld/8,
                                                                 conditional(x > 5*Ld/8,
                                                                             conditional(x < 3*Ld/4,
                                                                                         self.sigma * (x - 3*Ld/4)**2,
                                                                                         0.0),
                                                                             (x-Ld/2)**2+Ld**2/32),
                                                                 (x-Ld/4)**2),
                                                     0.0))
        elif self.Xi_family == 'proper_peak':
            if self.num_Xis > 1:
                raise NotImplementedError('Quadratic Xi not yet implemented for more than one Xi')
            else:
                self.Xi_functions.append(self.sigma * 0.5*2/(exp(x-Ld/2)+exp(-x+Ld/2)))

        elif self.Xi_family == 'constant':
            if self.num_Xis > 1:
                raise NotImplementedError('Constant Xi not yet implemented for more than one Xi')
            else:
                self.Xi_functions.append(self.sigma * (sin(0*pi*x/Ld)+1))


        else:
            raise NotImplementedError('Xi_family %s not implemented' % self.Xi_family)

        # make lists of functions for xi_x, xi_xx and xi_xxx
        if self.scheme in ['hydrodynamic', 'LASCH_hydrodynamic']:
            self.dXi_x = prognostic_variables.dXi_x
            self.dXi_xx = prognostic_variables.dXi_xx

            self.Xi_x_functions = []
            self.Xi_xx_functions = []

            for Xi_expr in self.Xi_functions:
                Xi_x_function = Function(self.dXi_x.function_space())
                Xi_xx_function = Function(self.dXi_xx.function_space())

                phi_x = TestFunction(self.dXi_x.function_space())
                phi_xx = TestFunction(self.dXi_xx.function_space())

                Xi_x_eqn = phi_x * Xi_x_function * dx + phi_x.dx(0) * Xi_expr * dx
                Xi_xx_eqn = phi_xx * Xi_xx_function * dx + phi_xx.dx(0) * Xi_x_function * dx

                Xi_x_problem = NonlinearVariationalProblem(Xi_x_eqn, Xi_x_function)
                Xi_xx_problem = NonlinearVariationalProblem(Xi_xx_eqn, Xi_xx_function)

                Xi_x_solver = NonlinearVariationalSolver(Xi_x_problem)
                Xi_xx_solver = NonlinearVariationalSolver(Xi_xx_problem)

                # for some reason these solvers don't work for constant Xi functions
                # so just manually make the derivatives be zero
                if self.Xi_family == 'constant':
                    Xi_x_function.interpolate(0.0*x)
                    Xi_xx_function.interpolate(0.0*x)
                else:
                    Xi_x_solver.solve()
                    Xi_xx_solver.solve()

                self.Xi_x_functions.append(Xi_x_function)
                self.Xi_xx_functions.append(Xi_xx_function)

        # now make a master xi
        Xi_expr = 0.0*x

        for dW, Xi_function, pure_xi, pure_xi_x, pure_xi_xx, pure_xi_xxx, pure_xi_xxxx in zip(self.dWs, self.Xi_functions, self.pure_xi_list, self.pure_xi_x_list, self.pure_xi_xx_list, self.pure_xi_xxx_list, self.pure_xi_xxxx_list):
            Xi_expr += dW * Xi_function
            if self.scheme in ['upwind', 'LASCH']:
                pure_xi.interpolate(as_vector([Xi_function]))
                pure_xi_x.project(as_vector([Xi_function.dx(0)]))

                CG1 = FunctionSpace(mesh, "CG", 1)
                psi =  TestFunction(CG1)
                xixx_scalar = Function(CG1)
                xixx_eqn = psi * xixx_scalar * dx + psi.dx(0) * Xi_function.dx(0) * dx
                prob = NonlinearVariationalProblem(xixx_eqn, xixx_scalar)
                solver = NonlinearVariationalSolver(prob)
                solver.solve()
                pure_xi_xx.interpolate(as_vector([xixx_scalar]))

            else:
                pure_xi.interpolate(Xi_function)

                # I guess we can't take the gradient of constants
                if self.Xi_family != 'constant':
                    pure_xi_x.project(Xi_function.dx(0))
                    pure_xi_xx.project(pure_xi_x.dx(0))
                    pure_xi_xxx.project(pure_xi_xx.dx(0))
                    pure_xi_xxxx.project(pure_xi_xxx.dx(0))

        if self.scheme in ['upwind', 'LASCH']:
            self.dXi_interpolator = Interpolator(as_vector([Xi_expr]), self.dXi)
        else:
            self.dXi_interpolator = Interpolator(Xi_expr, self.dXi)

        if self.scheme in ['hydrodynamic', 'LASCH_hydrodynamic']:

            # initialise blank expressions
            Xi_x_expr = 0.0*x
            Xi_xx_expr = 0.0*x

            # make full expressions by adding all dW * Xi_xs
            for dW, Xi_x_function, Xi_xx_function in zip(self.dWs, self.Xi_x_functions, self.Xi_xx_functions):
                Xi_x_expr += dW * Xi_x_function
                Xi_xx_expr += dW * Xi_xx_function

            self.dXi_x_interpolator = Interpolator(Xi_x_expr, self.dXi_x)
            self.dXi_xx_interpolator = Interpolator(Xi_xx_expr, self.dXi_xx)


    def update(self, t):
        """
        Updates the Xi function for the next time step.
        """
        if self.num_Xis > 0:
            # Try to calculate dW numbers separately

            # For nXi_updates > 1 then the ordering of calls to np.random.randn
            # needs to be equivalent to that for the corresponding smaller dt
            if self.nXi_updates > 1:
                self.dW_nums[:] = 0.0
                for j in range(self.nXi_updates):
                    for i in range(self.num_Xis):
                        self.dW_nums[i] += np.random.randn() * np.sqrt(self.dt/self.nXi_updates)

            else:
                for i in range(self.num_Xis):
                    if self.smooth_t is not None:
                        self.dW_nums[i] = self.smooth_t(t)
                    elif self.fixed_dW is not None:
                        self.dW_nums[i] = self.fixed_dW * np.sqrt(self.dt)
                    else:
                        self.dW_nums[i] = np.random.randn() * np.sqrt(self.dt)

            # This is to ensure we stick close to what the original code did
            [dw.assign(dw_num) for dw, dw_num in zip(self.dWs, self.dW_nums)]
            self.dXi_interpolator.interpolate()

            if self.scheme in ['hydrodynamic', 'LASCH_hydrodynamic']:
                self.dXi_x_interpolator.interpolate()
                self.dXi_xx_interpolator.interpolate()
        else:
            pass
