from firedrake import (Interpolator, Constant, as_vector, sin,
                       cos, exp, FunctionSpace, VectorFunctionSpace,
                       pi, SpatialCoordinate)
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
        self.dWs = [Constant(0.0) for dw in range(self.num_Xis)]
        self.Xi_functions = []

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
                if (n+1) % 2 == 0:
                    self.Xi_functions.append(sin(2*(n+1)*pi*x/Ld))
                else:
                    self.Xi_functions.append(cos(2*(n+1)*pi*x/Ld))

        elif self.Xi_family == 'gaussians':
            for n in range(self.num_Xis):
                    self.Xi_functions.append(0.5*self.num_Xis*exp(-((x-Ld*(n+1)/(self.num_Xis +1.0))/2.)**2))

        else:
            raise NotImplementedError('Xi_family %s not implemented' % self.Xi_family)

        Xi_expr = 0.0*x

        for dW, Xi_function in zip(self.dWs, self.Xi_functions):
            Xi_expr += dW * Xi_function

        if simulation_parameters['scheme'][-1] == 'upwind':
            self.Xi_interpolator = Interpolator(as_vector([Xi_expr]), self.Xi)
        else:
            self.Xi_interpolator = Interpolator(Xi_expr, self.Xi)



    def update(self):
        """
        Updates the Xi function for the next time step.
        """
        if self.num_Xis > 0:
            [dw.assign(self.sigma*np.random.randn()) for dw in self.dWs]
            self.Xi_interpolator.interpolate()
        else:
            pass
