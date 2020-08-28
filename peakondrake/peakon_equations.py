from firedrake import (SpatialCoordinate, Function)
import numpy as np

class PeakonEquations(object):
    """
    An object setting the stochastic (ordinary) differential equations
    describing the motion of a peakon in the stochastic Camassa-Holm equation.

    :arg prognostic_variables: a PrognosticVariables object.
    :arg simulation_parameters: a dictionary storing the simulation parameters.
    """

    def __init__(self, prognostic_variables, simulation_parameters, peakon_speed, method='ito_euler'):

        self.dWs = prognostic_variables.dW_nums
        self.pure_xi_list = prognostic_variables.pure_xi_list
        self.pure_xi_x_list = prognostic_variables.pure_xi_x_list
        self.pure_xi_xx_list = prognostic_variables.pure_xi_xx_list
        self.u = prognostic_variables.u
        mesh = simulation_parameters['mesh'][-1]
        x, = SpatialCoordinate(mesh)
        coords = Function(self.u.function_space()).project(x)
        self.dt = simulation_parameters['dt'][-1]
        self.Ld = simulation_parameters['Ld'][-1]
        self.method = method

        p0 = np.max(self.u.dat.data[:])
        q0 = coords.dat.data[np.argmax(self.u.dat.data[:])]

        self.p = peakon_speed if peakon_speed is not None else p0
        self.q = q0

    def update(self):

        dp = 0
        dq = self.p*self.dt

        if self.method == 'ito_euler':

            for xi_field, xi_x_field, xi_xx_field, dW in zip(self.pure_xi_list, self.pure_xi_x_list, self.pure_xi_xx_list, self.dWs):
                xi = xi_field.at(self.q, tolerance=1e-6)
                xi_x = xi_x_field.at(self.q, tolerance=1e-6)
                xi_xx = xi_xx_field.at(self.q, tolerance=1e-6)

                dp += self.p/2*(xi_x*xi_x - xi*xi_xx)*self.dt - self.p*xi_x*dW
                dq += 0.5*xi*xi_x*self.dt + xi*dW

        elif self.method == 'milstein':

            for xi_field, xi_x_field, xi_xx_field, dW in zip(self.pure_xi_list, self.pure_xi_x_list, self.pure_xi_xx_list, self.dWs):
                xi = xi_field.at(self.q, tolerance=1e-6)
                xi_x = xi_x_field.at(self.q, tolerance=1e-6)
                xi_xx = xi_xx_field.at(self.q, tolerance=1e-6)

                dp += self.p/2*(xi_x*xi_x - xi*xi_xx)*dW**2 - self.p*xi_x*dW
                dq += 0.5*xi*xi_x*dW**2 + xi*dW

        self.p += dp
        self.q += dq

        if self.q < 0:
            self.q += self.Ld
        if self.q > self.Ld:
            self.q -= self.Ld
