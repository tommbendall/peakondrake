from firedrake import (dx, dS, FacetNormal, TestFunction, Function,
                       Constant, jump, NonlinearVariationalSolver,
                       NonlinearVariationalProblem, dot, as_vector, sqrt,
                       MixedFunctionSpace, TestFunctions, split)

class Equations(object):
    """
    An object setting up equations and solvers for
    the stochastic Camassa-Holm equation.

    :arg prognostic_variables: a PrognosticVariables object.
    :arg simulation_parameters: a dictionary storing the simulation parameters.
    """

    def __init__(self, prognostic_variables, simulation_parameters):

        mesh = simulation_parameters['mesh'][-1]
        self.scheme = simulation_parameters['scheme'][-1]
        self.timestepping = simulation_parameters['timestepping'][-1]
        alphasq = simulation_parameters['alphasq'][-1]
        c0 = simulation_parameters['c0'][-1]
        gamma = simulation_parameters['gamma'][-1]
        Dt = Constant(simulation_parameters['dt'][-1])
        self.solvers = []

        if alphasq.values()[0] > 0.0 and gamma.values()[0] == 0.0:
            self.setup = 'ch'
            if self.scheme == 'upwind' and self.timestepping == 'ssprk3':

                Vm = prognostic_variables.Vm
                Vu = prognostic_variables.Vu
                self.m = prognostic_variables.m
                self.u = prognostic_variables.u
                self.Xi = prognostic_variables.Xi
                self.m0 = Function(Vm).assign(self.m)

                # now make problem for the actual problem
                psi = TestFunction(Vm)
                self.m_trial = Function(Vm)
                self.dm = Function(Vm)  # introduce this as the advection operator for a single step

                us = Dt * self.u + sqrt(Dt) * self.Xi

                nhat = FacetNormal(mesh)
                un = 0.5*(dot(us, nhat) + abs(dot(us, nhat)))
                ones = Function(Vu).project(as_vector([Constant(1.)]))

                Lm = (psi * self.dm * dx
                      - psi.dx(0) * self.m_trial * dot(ones, us) * dx
                      + psi* self.m_trial * dot(ones, us.dx(0)) * dx
                      + jump(psi) * (un('+')*self.m_trial('+') - un('-')*self.m_trial('-')) * dS)
                mprob = NonlinearVariationalProblem(Lm, self.dm)
                self.msolver = NonlinearVariationalSolver(mprob, solver_parameters={'ksp_type':'preonly',
                                                                                    'pc_type':'bjacobi',
                                                                                    'sub_pc_type':'ilu'})

                phi = TestFunction(Vu)
                Lu = (dot(phi, ones) * self.m * dx - dot(phi, self.u) * dx - alphasq * dot(self.u.dx(0), phi.dx(0)) * dx)
                uprob = NonlinearVariationalProblem(Lu, self.u)
                self.usolver = NonlinearVariationalSolver(uprob, solver_parameters={'ksp_type':'preonly',
                                                                                    'pc_type':'lu'})

            elif self.scheme == 'upwind' and self.timestepping == 'midpoint':
                Vm = prognostic_variables.Vm
                Vu = prognostic_variables.Vu
                self.m = prognostic_variables.m
                self.u = prognostic_variables.u
                self.Xi = prognostic_variables.Xi
                self.m0 = Function(Vm).assign(self.m)

                # now make problem for the actual problem
                psi = TestFunction(Vm)
                self.m_trial = Function(Vm)
                self.mh = (self.m0 + self.m_trial) / 2

                us = Dt * self.u + sqrt(Dt) * self.Xi

                nhat = FacetNormal(mesh)
                un = 0.5*(dot(us, nhat) + abs(dot(us, nhat)))
                ones = Function(Vu).project(as_vector([Constant(1.)]))

                Lm = (psi * self.m_trial * dx - psi * self.m0 * dx
                      - psi.dx(0) * self.mh * dot(ones, us) * dx
                      + psi* self.mh * dot(ones, us.dx(0)) * dx
                      + jump(psi) * (un('+')*self.mh('+') - un('-')*self.mh('-')) * dS)
                mprob = NonlinearVariationalProblem(Lm, self.m_trial)
                self.msolver = NonlinearVariationalSolver(mprob, solver_parameters={'ksp_type':'preonly',
                                                                                    'pc_type':'bjacobi',
                                                                                    'sub_pc_type':'ilu'})

                phi = TestFunction(Vu)
                Lu = (dot(phi, ones) * self.m * dx - dot(phi, self.u) * dx - alphasq * dot(self.u.dx(0), phi.dx(0)) * dx)
                uprob = NonlinearVariationalProblem(Lu, self.u)
                self.usolver = NonlinearVariationalSolver(uprob, solver_parameters={'ksp_type':'preonly',
                                                                                    'pc_type':'lu'})

            elif self.scheme == 'conforming' and self.timestepping == 'midpoint':
                Vm = prognostic_variables.Vm
                Vu = prognostic_variables.Vu

                self.u = prognostic_variables.u
                self.m = prognostic_variables.m
                self.Xi = prognostic_variables.Xi
                self.u0 = prognostic_variables.u0.assign(self.u)

                zeta = TestFunction(Vm)
                m_eqn = zeta * self.m * dx - zeta * self.u0 * dx - alphasq * zeta.dx(0) * self.u0.dx(0) * dx
                m_prob = NonlinearVariationalProblem(m_eqn, self.m)
                m_solver = NonlinearVariationalSolver(m_prob)
                m_solver.solve()

                W = MixedFunctionSpace((Vu, Vm))
                psi, phi = TestFunctions(W)

                w1 = Function(W)
                self.u1, self.m1 = split(w1)
                uh = (self.u1 + self.u) / 2
                mh = (self.m1 + self.m) / 2
                us = Dt * uh + sqrt(Dt) * self.Xi

                Lu = (psi * (self.m1 - self.m) * dx
                      - psi.dx(0) * mh * us * dx
                      + psi * mh * us.dx(0) * dx
                      - phi * self.m1 * dx
                      + phi * self.u1 * dx
                      + alphasq * phi.dx(0) * self.u1.dx(0) * dx)

                self.u1, self.m1 = w1.split()

                uprob = NonlinearVariationalProblem(Lu, w1)
                self.usolver = NonlinearVariationalSolver(uprob,
                                                          solver_parameters={'mat_type': 'aij',
                                                                             'ksp_type': 'preonly',
                                                                             'pc_type': 'lu'})


            elif self.scheme == 'hydrodynamic' and self.timestepping == 'midpoint':
                Vu = prognostic_variables.Vu

                self.u = prognostic_variables.u

                W = MixedFunctionSpace((Vu,)*3)
                psi, phi, zeta = TestFunctions(W)

                w1 = Function(W)
                self.u1, dFh, dGh = split(w1)

                uh  = (self.u1 + self.u) / 2
                dXi = sqrt(Dt) * prognostic_variables.Xi
                dXi_x = sqrt(Dt) * prognostic_variables.Xi_x
                dXi_xx = sqrt(Dt) * prognostic_variables.Xi_xx
                dvh = Dt * uh + dXi

                Lu = (psi * (self.u1 - self.u) * dx
                      + psi * uh.dx(0) * dvh * dx
                      - psi.dx(0) * dFh * dx
                      + psi * dGh * dx
                      + phi * dFh * dx + alphasq * phi.dx(0) * dFh.dx(0) * dx
                      - phi * uh * uh * Dt * dx - 0.5 * alphasq * phi * uh.dx(0) * uh.dx(0) * Dt * dx
                      + zeta * dGh * dx + alphasq * zeta.dx(0) * dGh.dx(0) * dx
                      - 2 * zeta * uh * dXi_x * dx - alphasq * zeta * uh.dx(0) * dXi_xx * dx)

                self.u1, dFh, dGh = w1.split()

                uprob = NonlinearVariationalProblem(Lu, w1)
                self.usolver = NonlinearVariationalSolver(uprob,
                                                          solver_parameters={'mat_type': 'aij',
                                                                             'ksp_type': 'preonly',
                                                                             'pc_type': 'lu'})

            elif self.scheme == 'test' and self.timestepping == 'midpoint':
                self.u = prognostic_variables.u
                Vu = prognostic_variables.Vu
                psi = TestFunction(Vu)
                self.u1 = Function(Vu)
                uh = (self.u1 + self.u) / 2
                dvh = Dt * uh + sqrt(Dt) * prognostic_variables.Xi

                eqn = (psi * (self.u1 - self.u) * dx - psi * uh * dvh.dx(0) * dx)
                prob = NonlinearVariationalProblem(eqn, self.u1)
                self.usolver = NonlinearVariationalSolver(prob,
                                                          solver_parameters={'mat_type': 'aij',
                                                                             'ksp_type': 'preonly',
                                                                             'pc_type': 'lu'})

            else:
                raise ValueError('Scheme %s and timestepping %s either not compatible or not recognised.' % (self.scheme, self.timestepping))

        elif alphasq.values()[0] == 0.0 and gamma.values()[0] > 0.0:
            self.setup = 'kdv'
            if self.scheme == 'upwind' and self.timestepping == 'ssprk3':
                raise NotImplementedError('Scheme %s and timestepping %s not yet implemented.' % (self.scheme, self.timestepping))

            elif self.scheme == 'upwind' and self.timestepping == 'midpoint':
                raise NotImplementedError('Scheme %s and timestepping %s not yet implemented.' % (self.scheme, self.timestepping))

            elif self.scheme == 'conforming' and self.timestepping == 'midpoint':
                Vf = prognostic_variables.Vf
                Vu = prognostic_variables.Vu

                self.u = prognostic_variables.u
                self.Xi = prognostic_variables.Xi
                self.u0 = prognostic_variables.u0

                W = MixedFunctionSpace((Vu, Vf))
                psi, phi = TestFunctions(W)

                w1 = Function(W)
                self.u1, Fh = split(w1)
                uh = (self.u1 + self.u) / 2
                us = Dt * uh + sqrt(Dt) * self.Xi

                Lu = (psi * (self.u1 - self.u) * dx
                      - 6 * psi.dx(0) * uh * us * dx
                      + 6 * psi * uh * us.dx(0) * dx
                      - gamma * psi.dx(0) * Fh * dx
                      - phi * Fh * dx - phi.dx(0) * us.dx(0) * dx)

                self.u1, Fh = w1.split()

                uprob = NonlinearVariationalProblem(Lu, w1)
                self.usolver = NonlinearVariationalSolver(uprob,
                                                          solver_parameters={'mat_type': 'aij',
                                                                             'ksp_type': 'preonly',
                                                                             'pc_type': 'lu'})


            elif self.scheme == 'hydrodynamic' and self.timestepping == 'midpoint':
                raise NotImplementedError('Scheme %s and timestepping %s not yet implemented.' % (self.scheme, self.timestepping))

            else:
                raise ValueError('Scheme %s and timestepping %s either not compatible or not recognised.' % (self.scheme, self.timestepping))

        else:
            raise NotImplementedError('Schemes for your values of alpha squared %.3f and gamma %.3f are not yet implemented.' % (alphasq, gamma))


    def solve(self):

        if self.scheme == 'upwind' and self.timestepping == 'ssprk3':

            # do three step RK method for m
            self.m_trial.assign(self.m0)
            self.msolver.solve()
            self.m_trial.assign(self.m0 + self.dm)
            self.msolver.solve()
            self.m_trial.assign(3./4*self.m0 + 1./4*(self.m_trial + self.dm))
            self.msolver.solve()
            self.m.assign(1./3*self.m0 + 2./3*(self.m_trial + self.dm))
            self.m0.assign(self.m)

            # now solve inverse problem for u
            self.usolver.solve()

        elif self.scheme == 'upwind' and self.timestepping == 'midpoint':
            self.msolver.solve()
            self.m.assign(self.m_trial)
            self.usolver.solve()
            self.m0.assign(self.m)

        elif self.scheme == 'conforming' and self.timestepping == 'midpoint':
            if self.setup == 'kdv':
                self.usolver.solve()
                self.u0.assign(self.u)
                self.u.assign(self.u1)
            elif self.setup == 'ch':
                self.usolver.solve()
                self.u.assign(self.u1)
                self.m.assign(self.m1)

        elif self.scheme == 'hydrodynamic' and self.timestepping == 'midpoint':
            self.usolver.solve()
            self.u.assign(self.u1)

        elif self.scheme == 'test' and self.timestepping == 'midpoint':
            self.usolver.solve()
            self.u.assign(self.u1)

        else:
            raise ValueError('Scheme %s and timestepping %s either not compatible or not recognised.' % (self.scheme, self.timestepping))
