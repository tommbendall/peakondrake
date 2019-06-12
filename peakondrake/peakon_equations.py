from firedrake import (SpatialCoordinate, Function)
import numpy as np

class PeakonEquations(object):
    """
    An object setting the stochastic (ordinary) differential equations
    describing the motion of a peakon in the stochastic Camassa-Holm equation.

    :arg prognostic_variables: a PrognosticVariables object.
    :arg simulation_parameters: a dictionary storing the simulation parameters.
    """

    def __init__(self, prognostic_variables, simulation_parameters):

        self.Xi = prognostic_variables.Xi
        self.Xi_x = prognostic_variables.Xi_x
        self.pure_xi_list = prognostic_variables.pure_xi_list
        self.pure_xi_x_list = prognostic_variables.pure_xi_x_list
        self.u = prognostic_variables.u
        mesh = simulation_parameters['mesh'][-1]
        x, = SpatialCoordinate(mesh)
        coords = Function(self.u.function_space()).project(x)
        self.dt = simulation_parameters['dt'][-1]
        self.Ld = simulation_parameters['Ld'][-1]

        p0 = np.max(self.u.dat.data[:])
        q0 = coords.dat.data[np.argmax(self.u.dat.data[:])]

        self.p = p0
        self.q = q0

    def update(self):

        dp = - self.dt ** 0.5 * self.p * self.Xi_x.at(self.q, tolerance=1e-6)
        dq = self.p * self.dt + self.dt ** 0.5 * self.Xi.at(self.q, tolerance=1e-6)
        for pure_xi, pure_xi_x in zip(self.pure_xi_list, self.pure_xi_x_list):
            dp += 0.5 * self.p * pure_xi_x.at(self.q, tolerance=1e-6) ** 2 * self.dt
            dq += 0.5 * pure_xi.at(self.q, tolerance=1e-6) * pure_xi_x.at(self.q, tolerance=1e-6) * self.dt

        self.p += dp
        self.q += dq

        if self.q < 0:
            self.q += self.Ld
        if self.q > self.Ld:
            self.q -= self.Ld
