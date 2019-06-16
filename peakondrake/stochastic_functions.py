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

        self.num_Xis = simulation_parameters['num_Xis'][-1]
        self.Xi_family = simulation_parameters['Xi_family'][-1]
        self.Xi = prognostic_variables.Xi
        self.dWs = [Constant(0.0) for dw in range(self.num_Xis)]
        self.Xi_functions = []
        self.nXi_updates = simulation_parameters['nXi_updates'][-1]
        self.smooth_t = simulation_parameters['smooth_t'][-1]

        if self.smooth_t is not None and self.nXi_updates > 1:
            raise ValueError('Prescribing forcing and including multiple Xi updates are not compatible.')

        if self.smooth_t is not None:
            print('WARNING: Remember to change sigma to sigma * sqrt(dt) with the prescribed forcing option.')

        if self.nXi_updates > 1:
            print('WARNING: Remember to change sigma to sigma / sqrt(nXi_updates) with nXi_updates.')


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
        for xi in range(self.num_Xis):
            self.pure_xi_list.append(Function(self.Xi.function_space()))
            self.pure_xi_x_list.append(Function(self.Xi.function_space()))


        if self.Xi_family == 'sines':
            for n in range(self.num_Xis):
                if (n+1) % 2 == 1:
                    self.Xi_functions.append(sin(2*(n+1)*pi*x/Ld))
                else:
                    self.Xi_functions.append(cos(2*(n+1)*pi*x/Ld))

        elif self.Xi_family == 'gaussians':
            for n in range(self.num_Xis):
                self.Xi_functions.append(0.5*self.num_Xis*exp(-((x-Ld*(n+1)/(self.num_Xis +1.0))/2.)**2))

        elif self.Xi_family == 'quadratic':
            if self.num_Xis > 1:
                raise NotImplementedError('Quadratic Xi not yet implemented for more than one Xi')
            else:
                self.Xi_functions.append(32/(Ld*Ld)*conditional(x > Ld/4,
                                                     conditional(x > 3*Ld/8,
                                                                 conditional(x > 5*Ld/8,
                                                                             conditional(x < 3*Ld/4,
                                                                                         (x - 3*Ld/4)**2,
                                                                                         0.0),
                                                                             (x-Ld/2)**2+Ld**2/32),
                                                                 (x-Ld/4)**2),
                                                     0.0))

        else:
            raise NotImplementedError('Xi_family %s not implemented' % self.Xi_family)

        # make lists of functions for xi_x, xi_xx and xi_xxx
        if self.scheme == 'hydrodynamic':
            self.Xi_x = prognostic_variables.Xi_x
            self.Xi_xx = prognostic_variables.Xi_xx

            self.Xi_x_functions = []
            self.Xi_xx_functions = []

            for Xi_expr in self.Xi_functions:
                Xi_x_function = Function(self.Xi_x.function_space())
                Xi_xx_function = Function(self.Xi_xx.function_space())

                phi_x = TestFunction(self.Xi_x.function_space())
                phi_xx = TestFunction(self.Xi_xx.function_space())

                Xi_x_eqn = phi_x * Xi_x_function * dx + phi_x.dx(0) * Xi_expr * dx
                Xi_xx_eqn = phi_xx * Xi_xx_function * dx + phi_xx.dx(0) * Xi_x_function * dx

                Xi_x_problem = NonlinearVariationalProblem(Xi_x_eqn, Xi_x_function)
                Xi_xx_problem = NonlinearVariationalProblem(Xi_xx_eqn, Xi_xx_function)

                Xi_x_solver = NonlinearVariationalSolver(Xi_x_problem)
                Xi_xx_solver = NonlinearVariationalSolver(Xi_xx_problem)

                Xi_x_solver.solve()
                Xi_xx_solver.solve()

                self.Xi_x_functions.append(Xi_x_function)
                self.Xi_xx_functions.append(Xi_xx_function)

        # now make a master xi
        Xi_expr = 0.0*x

        for dW, Xi_function, pure_xi, pure_xi_x in zip(self.dWs, self.Xi_functions, self.pure_xi_list, self.pure_xi_x_list):
            Xi_expr += dW * Xi_function
            if self.scheme == 'upwind':
                pure_xi.interpolate(as_vector([self.sigma*Xi_function]))
                pure_xi_x.project(as_vector([self.sigma*Xi_function.dx(0)]))
            else:
                pure_xi.interpolate(self.sigma*Xi_function)
                pure_xi_x.project(self.sigma*Xi_function.dx(0))

        if self.scheme == 'upwind':
            self.Xi_interpolator = Interpolator(as_vector([Xi_expr]), self.Xi)
        else:
            self.Xi_interpolator = Interpolator(Xi_expr, self.Xi)

        if self.scheme == 'hydrodynamic':

            # initialise blank expressions
            Xi_x_expr = 0.0*x
            Xi_xx_expr = 0.0*x

            # make full expressions by adding all dW * Xi_xs
            for dW, Xi_x_function, Xi_xx_function in zip(self.dWs, self.Xi_x_functions, self.Xi_xx_functions):
                Xi_x_expr += dW * Xi_x_function
                Xi_xx_expr += dW * Xi_xx_function

            self.Xi_x_interpolator = Interpolator(Xi_x_expr, self.Xi_x)
            self.Xi_xx_interpolator = Interpolator(Xi_xx_expr, self.Xi_xx)


    def update(self, t):
        """
        Updates the Xi function for the next time step.
        """
        if self.num_Xis > 0:
            if self.smooth_t is not None:
                [dw.assign(self.sigma*self.smooth_t(t)) for dw in self.dWs]
            else:
                [dw.assign(self.sigma*np.random.randn()) for dw in self.dWs]
                if self.nXi_updates > 1:
                    for i in range(self.nXi_updates-1):
                        [dw.assign(dw + self.sigma*np.random.randn()) for dw in self.dWs]
            self.Xi_interpolator.interpolate()

            if self.scheme == 'hydrodynamic':
                self.Xi_x_interpolator.interpolate()
                self.Xi_xx_interpolator.interpolate()
        else:
            pass
