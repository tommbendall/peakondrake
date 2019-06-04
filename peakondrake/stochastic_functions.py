from firedrake import (Interpolator, Constant, as_vector, sin,
                       cos, exp, FunctionSpace, VectorFunctionSpace,
                       pi, SpatialCoordinate, Function, conditional)
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

        self.num_Xis = simulation_parameters['num_Xis'][-1]
        self.Xi_family = simulation_parameters['Xi_family'][-1]
        self.Xi = prognostic_variables.Xi
        self.pure_xis = prognostic_variables.pure_xis
        for xi in range(self.num_Xis):
            self.pure_xis.append(Function(self.Xi.function_space()))
        self.dWs = [Constant(0.0) for dw in range(self.num_Xis)]
        self.Xi_functions = []
        self.sigma_kick = simulation_parameters['sigma_kick'][-1]
        self.t_kick = simulation_parameters['t_kick'][-1]

        seed = simulation_parameters['seed'][-1]
        np.random.seed(seed)

        # make sure sigma is a Constant
        if self.num_Xis != 0:
            if isinstance(simulation_parameters['sigma'][-1], Constant):
                self.sigma = simulation_parameters['sigma'][-1] / self.num_Xis
            else:
                self.sigma = Constant(simulation_parameters['sigma'][-1] / self.num_Xis)
        else:
            self.sigma = Constant(0.0)

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

        Xi_expr = 0.0*x

        for dW, Xi_function, pure_xi in zip(self.dWs, self.Xi_functions, self.pure_xis):
            Xi_expr += dW * Xi_function
            if simulation_parameters['scheme'][-1] == 'upwind':
                pure_xi.interpolate(as_vector([Xi_function]))
            else:
                pure_xi.interpolate(Xi_function)

        if simulation_parameters['scheme'][-1] == 'upwind':
            self.Xi_interpolator = Interpolator(as_vector([Xi_expr]), self.Xi)
        else:
            self.Xi_interpolator = Interpolator(Xi_expr, self.Xi)



    def update(self, t):
        """
        Updates the Xi function for the next time step.
        """
        if self.num_Xis > 0:
            if len(self.t_kick) > 0 and np.min(self.t_kick) < t < np.max(self.t_kick):
                [dw.assign(self.sigma_kick*np.random.randn()) for dw in self.dWs]
            else:
                [dw.assign(self.sigma*np.random.randn()) for dw in self.dWs]
            self.Xi_interpolator.interpolate()
        else:
            pass
