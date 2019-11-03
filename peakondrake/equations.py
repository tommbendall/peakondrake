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

            elif self.scheme == 'LASCH' and self.timestepping == 'ssprk3':
                Vm = prognostic_variables.Vm
                Vu = prognostic_variables.Vu
                VL = prognostic_variables.VL
                self.u = prognostic_variables.u
                self.m = prognostic_variables.m
                self.Em = prognostic_variables.Em
                self.Eu = prognostic_variables.Eu
                self.Xi = prognostic_variables.Xi
                self.m0 = Function(Vm).assign(self.m)
                self.Em0 = Function(Vm).assign(self.Em)
                pure_xis = prognostic_variables.pure_xi_list
                pure_xixs = prognostic_variables.pure_xi_x_list
                pure_xixxs = prognostic_variables.pure_xi_xx_list

                Lxis = []
                L2xis = []
                self.Lxisolvers = []
                self.L2xisolvers = []
                for xi in pure_xis:
                    Lxis.append(Function(VL))
                    L2xis.append(Function(VL))

                # now make problem for the actual problem
                psi = TestFunction(Vm)
                self.m_trial = Function(Vm)
                self.dm = Function(Vm)  # introduce this as the advection operator for a single step
                self.Em_trial = Function(Vm)
                self.dEm = Function(Vm)

                us = Dt * self.Eu + sqrt(Dt) * self.Xi

                nhat = FacetNormal(mesh)
                usn = 0.5*(dot(us, nhat) + abs(dot(us, nhat)))
                ones = Function(Vu).project(as_vector([Constant(1.)]))

                Lm = (psi * self.dm * dx
                      - psi.dx(0) * self.m_trial * dot(ones, us) * dx
                      + psi* self.m_trial * dot(ones, us.dx(0)) * dx
                      + jump(psi) * (usn('+')*self.m_trial('+') - usn('-')*self.m_trial('-')) * dS)
                mprob = NonlinearVariationalProblem(Lm, self.dm)
                self.msolver = NonlinearVariationalSolver(mprob, solver_parameters={'ksp_type':'preonly',
                                                                                    'pc_type':'bjacobi',
                                                                                    'sub_pc_type':'ilu'})
                Eun = 0.5*(dot(self.Eu, nhat) + abs(dot(self.Eu, nhat)))
                nu = TestFunction(Vm)
                Emeqn = (nu * self.dEm * dx
                         - Dt * nu.dx(0) * self.Em_trial * dot(ones, self.Eu) * dx
                         + Dt * nu * self.Em_trial * dot(ones, self.Eu.dx(0)) * dx
                         + Dt * jump(nu) * (Eun('+')*self.Em_trial('+') - Eun('-')*self.Em_trial('-')) * dS)
                for xi, xix, xixx, Lxi, L2xi in zip(pure_xis, pure_xixs, pure_xixxs, Lxis, L2xis):
                    # xin = 0.5*(dot(dot(ones, xi)*xi, nhat) + abs(dot(dot(ones, xi)*xi, nhat)))
                    # xixn = 0.5*(dot(dot(ones, xi)*xix, nhat) + abs(dot(dot(ones, xi)*xix, nhat)))
                    # Emeqn -= 0.5 * Dt * (nu * self.Em_trial * (2 * dot(ones, xi) * dot(ones, xixx) + 4 * dot(ones, xix) * dot(ones, xix)) * dx
                    #                      - 5 * (nu.dx(0) * self.Em_trial * dot(ones, xix) * dot(ones, xi) * dx
                    #                             + nu * self.Em_trial * dot(ones, xix.dx(0)) * dot(ones, xi) * dx
                    #                             + nu * self.Em_trial * dot(ones, xix) * dot(ones, xi.dx(0)) * dx
                    #                             - jump(nu) * (xixn('+')*self.Em_trial('+') - xixn('-')*self.Em_trial('-')) * dS)
                    #                      - (nu.dx(0) * self.Em_trial.dx(0) * dot(ones, xi) * dot(ones, xi) * dx
                    #                         + 2 * nu * self.Em_trial.dx(0) * dot(ones, xi) * dot(ones, xix) * dx
                    #                         - jump(nu) * (xin('+')*self.Em_trial.dx(0)('+') - xin('-')*self.Em_trial.dx(0)('-')) * dS))
                    # L
                    # solve for Lk
                    chi = TestFunction(VL)
                    Leqn = (chi * Lxi * dx
                            + chi * self.Em_trial * dot(ones, xi) * dx
                            - chi * self.Em_trial * dot(ones, xi.dx(0)) * dx)
                    Lprob = NonlinearVariationalProblem(Leqn, Lxi)
                    Lsolver = NonlinearVariationalSolver(Lprob)
                    self.Lxisolvers.append(Lsolver)
                    # then solve for LLk
                    iota = TestFunction(VL)
                    L2eqn = (iota * L2xi * dx
                             + iota * Lxi * dot(ones, xi) * dx
                             - iota * Lxi * dot(ones, xi.dx(0)) * dx)
                    L2prob = NonlinearVariationalProblem(L2eqn, L2xi)
                    L2solver = NonlinearVariationalSolver(L2prob)
                    self.L2xisolvers.append(L2solver)

                    Emeqn -= 0.5 * Dt * nu * L2xi * dx
                Emprob = NonlinearVariationalProblem(Emeqn, self.dEm)
                self.Emsolver = NonlinearVariationalSolver(Emprob, solver_parameters={'ksp_type':'preonly',
                                                                                   'pc_type':'bjacobi',
                                                                                   'sub_pc_type':'ilu'})

                phi = TestFunction(Vu)
                Eueqn = (dot(phi, ones) * self.Em * dx - dot(phi, self.Eu) * dx - alphasq * dot(self.Eu.dx(0), phi.dx(0)) * dx)
                Euprob = NonlinearVariationalProblem(Eueqn, self.Eu)
                self.Eusolver = NonlinearVariationalSolver(Euprob, solver_parameters={'ksp_type':'preonly',
                                                                                    'pc_type':'lu'})

                zeta = TestFunction(Vu)
                ueqn = (dot(zeta, ones) * self.m * dx - dot(zeta, self.u) * dx - alphasq * dot(self.u.dx(0), zeta.dx(0)) * dx)
                uprob = NonlinearVariationalProblem(ueqn, self.u)
                self.usolver = NonlinearVariationalSolver(uprob, solver_parameters={'ksp_type':'preonly',
                                                                                    'pc_type':'lu'})

            elif self.scheme == 'LASCH' and self.timestepping == 'midpoint':
                Vm = prognostic_variables.Vm
                Vu = prognostic_variables.Vu
                self.u = prognostic_variables.u
                self.m = prognostic_variables.m
                self.Em = prognostic_variables.Em
                self.Eu = prognostic_variables.Eu
                self.Xi = prognostic_variables.Xi
                self.m0 = Function(Vm).assign(self.m)
                self.Em0 = Function(Vm).assign(self.Em)
                pure_xis = prognostic_variables.pure_xi_list
                pure_xixs = prognostic_variables.pure_xi_x_list
                pure_xixxs = prognostic_variables.pure_xi_xx_list

                # now make problem for the actual problem
                psi = TestFunction(Vm)
                self.m_trial = Function(Vm)
                self.mh = 0.5*(self.m_trial + self.m0)
                self.Em_trial = Function(Vm)
                self.Emh = 0.5*(self.Em_trial + self.Em0)

                us = Dt * self.Eu + sqrt(Dt) * self.Xi

                nhat = FacetNormal(mesh)
                usn = 0.5*(dot(us, nhat) + abs(dot(us, nhat)))
                ones = Function(Vu).project(as_vector([Constant(1.)]))

                Lm = (psi * self.m_trial * dx - psi * self.Em0 * dx
                      - psi.dx(0) * self.mh * dot(ones, us) * dx
                      + psi* self.mh * dot(ones, us.dx(0)) * dx
                      + jump(psi) * (usn('+')*self.mh('+') - usn('-')*self.mh('-')) * dS)
                mprob = NonlinearVariationalProblem(Lm, self.m_trial)
                self.msolver = NonlinearVariationalSolver(mprob, solver_parameters={'ksp_type':'preonly',
                                                                                    'pc_type':'bjacobi',
                                                                                    'sub_pc_type':'ilu'})
                Eun = 0.5*(dot(self.Eu, nhat) + abs(dot(self.Eu, nhat)))
                nu = TestFunction(Vm)
                Emeqn = (nu * self.Em_trial * dx - nu * self.Em0 * dx
                         - Dt * nu.dx(0) * self.Emh * dot(ones, self.Eu) * dx
                         + Dt * nu* self.Emh * dot(ones, self.Eu.dx(0)) * dx
                         + Dt * jump(nu) * (Eun('+')*self.Emh('+') - Eun('-')*self.Emh('-')) * dS)
                for xi, xix, xixx in zip(pure_xis, pure_xixs, pure_xixxs):
                    xin = 0.5*(dot(dot(ones, xi)*xi, nhat) + abs(dot(dot(ones, xi)*xi, nhat)))
                    xixn = 0.5*(dot(dot(ones, xi)*xix, nhat) + abs(dot(dot(ones, xi)*xix, nhat)))
                    Emeqn -= 0.5 * Dt * (nu * self.Emh * (2 * dot(ones, xi) * dot(ones, xixx) + 4 * dot(ones, xix) * dot(ones, xix)) * dx
                                         - 5 * (nu.dx(0) * self.Emh * dot(ones, xix) * dot(ones, xi) * dx
                                                + nu * self.Emh * dot(ones, xixx) * dot(ones, xi) * dx
                                                + nu * self.Emh * dot(ones, xix) * dot(ones, xix) * dx
                                                - jump(nu) * (xixn('+')*self.Emh('+') - xixn('-')*self.Emh('-')) * dS)
                                         - (nu.dx(0) * self.Emh.dx(0) * dot(ones, xi) * dot(ones, xi) * dx
                                            + 2 * nu * self.Emh.dx(0) * dot(ones, xi) * dot(ones, xix) * dx
                                            - jump(nu) * (xin('+')*self.Emh.dx(0)('+') - xin('-')*self.Emh.dx(0)('-')) * dS))

                Emprob = NonlinearVariationalProblem(Emeqn, self.Em_trial)
                self.Emsolver = NonlinearVariationalSolver(Emprob, solver_parameters={'ksp_type':'preonly',
                                                                                   'pc_type':'bjacobi',
                                                                                   'sub_pc_type':'ilu'})

                phi = TestFunction(Vu)
                Eueqn = (dot(phi, ones) * self.Em * dx - dot(phi, self.Eu) * dx - alphasq * dot(self.Eu.dx(0), phi.dx(0)) * dx)
                Euprob = NonlinearVariationalProblem(Eueqn, self.Eu)
                self.Eusolver = NonlinearVariationalSolver(Euprob, solver_parameters={'ksp_type':'preonly',
                                                                                    'pc_type':'lu'})

                zeta = TestFunction(Vu)
                ueqn = (dot(zeta, ones) * self.m * dx - dot(zeta, self.u) * dx - alphasq * dot(self.u.dx(0), zeta.dx(0)) * dx)
                uprob = NonlinearVariationalProblem(ueqn, self.u)
                self.usolver = NonlinearVariationalSolver(uprob, solver_parameters={'ksp_type':'preonly',
                                                                                    'pc_type':'lu'})

            elif self.scheme == 'LASCH_hydrodynamic' and self.timestepping == 'midpoint':
                Vu = prognostic_variables.Vu
                self.u = prognostic_variables.u
                self.Eu = prognostic_variables.Eu
                self.Xi = prognostic_variables.Xi
                pure_xis = prognostic_variables.pure_xi_list
                pure_xixs = prognostic_variables.pure_xi_x_list
                pure_xixxs = prognostic_variables.pure_xi_xx_list
                pure_xixxxs = prognostic_variables.pure_xi_xxx_list
                pure_xixxxxs = prognostic_variables.pure_xi_xxxx_list

                We = MixedFunctionSpace((Vu,)*3)

                Epsi, Ephi, Ezeta = TestFunctions(We)
                Ew1 = Function(We)
                self.Eu1, EFh, EGh = split(Ew1)

                Euh = (self.Eu1 + self.Eu) / 2

                Eueqn = (Epsi * (self.Eu1 - self.Eu) * dx
                         + Dt * Epsi * Euh * Euh.dx(0) * dx
                         - Dt * Epsi.dx(0) * EFh * dx
                         - Dt * Epsi * EGh * dx
                         + Ephi * EFh * dx + alphasq * Ephi.dx(0) * EFh.dx(0) * dx
                         - Ephi * 0.5 * alphasq * Euh.dx(0) * Euh.dx(0) * dx - Ephi * Euh * Euh * dx
                         + Ezeta * EGh * dx + alphasq * Ezeta.dx(0) * EGh.dx(0) * dx)

                for xi, xix, xixx, xixxx, xixxxx in zip(pure_xis, pure_xixs, pure_xixxs, pure_xixxxs, pure_xixxxxs):
                    Eueqn += (0.5 * Dt * Epsi.dx(0) * Euh.dx(0) * xi * xi * dx + Dt * Epsi.dx(0) * Euh * xi * xix * dx
                              + 1.5 * Dt * Epsi * Euh.dx(0) * xi * xix * dx + 0.5 * Dt * Epsi * Euh * (3 * xixx * xi + 2 * xix * xix) * dx
                              - 0.5 * Ezeta * Euh * (4 * xix * xix + 4 * xi * xixx - alphasq * xi * xixxxx - 2 * alphasq * xix * xixxx - alphasq * xixx * xixx) * dx
                              - 0.5 * Ezeta * Euh.dx(0) * (3 * xi * xix - alphasq * xixxx * xi + alphasq * xixx * xix) * dx)

                self.Eu1, EFh, EGh = Ew1.split()

                W = MixedFunctionSpace((Vu,)*4)

                psi, phi, zeta, Lambda = TestFunctions(W)
                w1 = Function(W)
                self.u1, dFh, dGh, H = split(w1)
                uh = (self.u1 + self.u) / 2
                dXi = sqrt(Dt) * prognostic_variables.Xi
                dXi_x = sqrt(Dt) * prognostic_variables.Xi_x
                dXi_xx = sqrt(Dt) * prognostic_variables.Xi_xx
                dv = Dt * (self.Eu + self.Eu1) / 2 + dXi
                DU = (self.Eu + self.Eu1) / 2 - uh

                ueqn = (psi * (self.u1 - self.u) * dx
                        + psi * uh.dx(0) * dv * dx
                        - psi.dx(0) * dFh * dx
                        + psi * dGh * dx
                        + phi * dFh * dx + alphasq * phi.dx(0) * dFh.dx(0) * dx
                        - phi * Dt * uh * uh * dx - alphasq / 2 * phi * Dt * uh.dx(0) * uh.dx(0) * dx
                        + zeta * dGh * dx + alphasq * zeta.dx(0) * dGh.dx(0) * dx
                        - 2 * zeta * uh * dXi.dx(0) * dx - alphasq * zeta * uh.dx(0) * dXi_xx * dx
                        - 2 * zeta * Dt * uh * DU.dx(0) * dx - alphasq * zeta * Dt * uh.dx(0) * H * dx
                        + Lambda * H * dx + Lambda.dx(0) * DU.dx(0) * dx)

                self.u1, dFh, dGh, H = w1.split()

                uprob = NonlinearVariationalProblem(ueqn, w1)
                Euprob = NonlinearVariationalProblem(Eueqn, Ew1)
                self.usolver = NonlinearVariationalSolver(uprob,
                                                          solver_parameters={'mat_type': 'aij',
                                                                             'ksp_type': 'preonly',
                                                                             'pc_type': 'lu'})
                self.Eusolver = NonlinearVariationalSolver(Euprob,
                                                           solver_parameters={'mat_type': 'aij',
                                                                             'ksp_type': 'preonly',
                                                                             'pc_type': 'lu'})

            elif self.scheme == 'LASCH_hydrodynamic_m' and self.timestepping == 'midpoint':
                Vu = prognostic_variables.Vu
                self.u = prognostic_variables.u
                self.Eu = prognostic_variables.Eu
                self.Xi = prognostic_variables.Xi
                pure_xis = prognostic_variables.pure_xi_list
                pure_xixs = prognostic_variables.pure_xi_x_list
                pure_xixxs = prognostic_variables.pure_xi_xx_list
                pure_xixxxs = prognostic_variables.pure_xi_xxx_list
                pure_xixxxxs = prognostic_variables.pure_xi_xxxx_list

                self.Emh = prognostic_variables.Em


                We = MixedFunctionSpace((Vu,)*3)

                Epsi, Ephi, Ezeta = TestFunctions(We)
                Ew1 = Function(We)
                self.Eu1, EFh, EGh = split(Ew1)

                Euh = (self.Eu1 + self.Eu) / 2

                Eueqn = (Epsi * (self.Eu1 - self.Eu) * dx
                         + Dt * Epsi * Euh * Euh.dx(0) * dx
                         - Dt * Epsi.dx(0) * EFh * dx
                         - Dt * Epsi * EGh * dx
                         + Ephi * EFh * dx + alphasq * Ephi.dx(0) * EFh.dx(0) * dx
                         - Ephi * 0.5 * alphasq * Euh.dx(0) * Euh.dx(0) * dx - Ephi * Euh * Euh * dx
                         + Ezeta * EGh * dx + alphasq * Ezeta.dx(0) * EGh.dx(0) * dx)

                for xi, xix, xixx, xixxx, xixxxx in zip(pure_xis, pure_xixs, pure_xixxs, pure_xixxxs, pure_xixxxxs):
                    Eueqn += (0.5 * Dt * Epsi.dx(0) * Euh.dx(0) * xi * xi * dx + Dt * Epsi.dx(0) * Euh * xi * xix * dx
                              + 1.5 * Dt * Epsi * Euh.dx(0) * xi * xix * dx + 0.5 * Dt * Epsi * Euh * (3 * xixx * xi + 2 * xix * xix) * dx
                              - 0.5 * Ezeta * Euh * (4 * xix * xix + 4 * xi * xixx - alphasq * xi * xixxxx - 2 * alphasq * xix * xixxx - alphasq * xixx * xixx) * dx
                              - 0.5 * Ezeta * Euh.dx(0) * (3 * xi * xix - alphasq * xixxx * xi + alphasq * xixx * xix) * dx)

                self.Eu1, EFh, EGh = Ew1.split()
                Euh = (self.Eu1 + self.Eu) / 2

                mu = TestFunction(Vu)

                # Emeqn = mu * self.Emh * dx - mu * Euh * dx + alphasq * mu.dx(0) * Euh.dx(0) * dx
                # Emprob = NonlinearVariationalProblem(Emeqn, self.Emh)
                Euh_xx = Function(Vu)
                Emeqn = mu * Euh_xx * dx + mu.dx(0) * Euh.dx(0) * dx
                Emprob = NonlinearVariationalProblem(Emeqn, Euh_xx)
                self.Emsolver = NonlinearVariationalSolver(Emprob)

                W = MixedFunctionSpace((Vu,)*2)

                psi, phi = TestFunctions(W)
                w1 = Function(W)
                self.u1, dFh = split(w1)
                uh = (self.u1 + self.u) / 2
                dXi = sqrt(Dt) * prognostic_variables.Xi
                dXi_x = sqrt(Dt) * prognostic_variables.Xi_x
                dXi_xx = sqrt(Dt) * prognostic_variables.Xi_xx
                dv = Dt * Euh + dXi

                ueqn = (psi * (self.u1 - self.u) * dx
                        + psi * uh.dx(0) * dv * dx
                        + psi * dFh * dx
                        + phi * dFh * dx + alphasq * phi.dx(0) * dFh.dx(0) * dx
                        - 2 * phi * Dt * uh * Euh.dx(0) * dx
                        - alphasq * phi * uh.dx(0) * Euh_xx * Dt * dx
                        - phi * (2 * uh * dXi_x + alphasq * uh.dx(0) * dXi_xx) * dx)

                self.u1, dFh = w1.split()

                uprob = NonlinearVariationalProblem(ueqn, w1)
                Euprob = NonlinearVariationalProblem(Eueqn, Ew1)
                self.usolver = NonlinearVariationalSolver(uprob,
                                                          solver_parameters={'mat_type': 'aij',
                                                                             'ksp_type': 'preonly',
                                                                             'pc_type': 'lu',
                                                                             'ksp_rtol': 1.0e-8,
                                                                             'ksp_atol': 1.0e-8})
                self.Eusolver = NonlinearVariationalSolver(Euprob,
                                                           solver_parameters={'mat_type': 'aij',
                                                                             'ksp_type': 'preonly',
                                                                             'pc_type': 'lu'})



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

            elif self.scheme == 'no_gradient' and self.timestepping == 'midpoint':
                # a version of the hydrodynamic form but without exploiting the gradient
                Vu = prognostic_variables.Vu

                self.u = prognostic_variables.u

                W = MixedFunctionSpace((Vu,)*3)
                psi, phi, zeta = TestFunctions(W)

                w1 = Function(W)
                self.u1, dFh, uh_xx = split(w1)

                uh  = (self.u1 + self.u) / 2
                dXi = sqrt(Dt) * prognostic_variables.Xi
                dXi_x = sqrt(Dt) * prognostic_variables.Xi_x
                dXi_xx = sqrt(Dt) * prognostic_variables.Xi_xx
                dvh = Dt * uh + dXi

                Lu = (psi * (self.u1 - self.u) * dx
                      + psi * uh.dx(0) * dvh * dx
                      + psi * dFh * dx
                      + phi * dFh * dx + alphasq * phi.dx(0) * dFh.dx(0) * dx
                      - 2 * phi * uh * uh.dx(0) * Dt * dx - alphasq * phi * uh.dx(0) * uh_xx * Dt * dx
                      - phi * (2 * uh * dXi_x + alphasq * uh.dx(0) * dXi_xx) * dx
                      + zeta * uh_xx * dx + zeta.dx(0) * uh.dx(0) * dx)

                self.u1, dFh, uh_xx = w1.split()

                uprob = NonlinearVariationalProblem(Lu, w1)
                self.usolver = NonlinearVariationalSolver(uprob,
                                                          solver_parameters={'mat_type': 'aij',
                                                                             'ksp_type': 'preonly',
                                                                             'pc_type': 'lu',
                                                                             'ksp_rtol': 1.0e-8,
                                                                             'ksp_atol': 1.0e-8})

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

        elif self.scheme == 'LASCH' and self.timestepping == 'ssprk3':
            # do three step RK method for m
            self.m_trial.assign(self.m0)
            self.msolver.solve()
            self.m_trial.assign(self.m0 + self.dm)
            self.msolver.solve()
            self.m_trial.assign(3./4*self.m0 + 1./4*(self.m_trial + self.dm))
            self.msolver.solve()
            self.m.assign(1./3*self.m0 + 2./3*(self.m_trial + self.dm))
            self.m0.assign(self.m)

            # do three step RK method for Em
            self.Em_trial.assign(self.Em0)
            for Lsolver, L2solver in zip(self.Lxisolvers, self.L2xisolvers):
                Lsolver.solve()
                L2solver.solve()
            self.Emsolver.solve()
            self.Em_trial.assign(self.Em0 + self.dEm)
            for Lsolver, L2solver in zip(self.Lxisolvers, self.L2xisolvers):
                Lsolver.solve()
                L2solver.solve()
            self.Emsolver.solve()
            self.Em_trial.assign(3./4*self.Em0 + 1./4*(self.Em_trial + self.dEm))
            for Lsolver, L2solver in zip(self.Lxisolvers, self.L2xisolvers):
                Lsolver.solve()
                L2solver.solve()
            self.Emsolver.solve()
            self.Em.assign(1./3*self.Em0 + 2./3*(self.Em_trial + self.dEm))
            self.Em0.assign(self.Em)

            # now solve inverse problem for u
            self.Eusolver.solve()
            self.usolver.solve()

        elif self.scheme == 'LASCH' and self.timestepping == 'midpoint':
            self.msolver.solve()
            self.m.assign(self.m_trial)
            self.Emsolver.solve()
            self.Em.assign(self.Em_trial)
            self.usolver.solve()
            self.Eusolver.solve()
            self.m0.assign(self.m)
            self.Em0.assign(self.Em)

        elif self.scheme == 'LASCH_hydrodynamic' and self.timestepping == 'midpoint':
            self.Eusolver.solve()
            self.usolver.solve()
            self.Eu.assign(self.Eu1)
            self.u.assign(self.u1)

        elif self.scheme == 'LASCH_hydrodynamic_m' and self.timestepping == 'midpoint':
            self.Eusolver.solve()
            self.Emsolver.solve()
            self.usolver.solve()
            self.Eu.assign(self.Eu1)
            self.u.assign(self.u1)

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

        elif self.scheme == 'no_gradient' and self.timestepping == 'midpoint':
            self.usolver.solve()
            self.u.assign(self.u1)

        elif self.scheme == 'test' and self.timestepping == 'midpoint':
            self.usolver.solve()
            self.u.assign(self.u1)

        else:
            raise ValueError('Scheme %s and timestepping %s either not compatible or not recognised.' % (self.scheme, self.timestepping))
