import numpy as np
from firedrake import (dx, conditional, exp, as_vector, dot, pi, sqrt,
                       Function, NonlinearVariationalSolver, cos,
                       TestFunction, NonlinearVariationalProblem,
                       SpatialCoordinate, Constant, FunctionSpace, cosh,
                       MixedFunctionSpace, TestFunctions, split)

def build_initial_conditions(prognostic_variables, simulation_parameters):

    """
    Initialises the prognostic variables based on the
    initial condition string.

    :arg prognostic_variables: a PrognosticVariables object.
    :arg simulation_parameters: a dictionary containing the simulation parameters.
    """

    mesh = simulation_parameters['mesh'][-1]
    ic = simulation_parameters['ic'][-1]
    alphasq = simulation_parameters['alphasq'][-1]
    c0 = simulation_parameters['c0'][-1]
    gamma = simulation_parameters['gamma'][-1]
    x, = SpatialCoordinate(mesh)
    Ld = simulation_parameters['Ld'][-1]
    deltax = Ld / simulation_parameters['resolution'][-1]
    epsilon = 1

    ic_dict = {'two_peaks': (0.2*2/(exp(x-403./15.*40./Ld) + exp(-x+403./15.*40./Ld))
                             + 0.5*2/(exp(x-203./15.*40./Ld)+exp(-x+203./15.*40./Ld))),
               'gaussian': 0.5*exp(-((x-10.)/2.)**2),
               'gaussian_narrow': 0.5*exp(-((x-10.)/1.)**2),
               'gaussian_wide': 0.5*exp(-((x-10.)/3.)**2),
               'peakon': conditional(x < Ld/2., exp((x-Ld/2)/sqrt(alphasq)), exp(-(x-Ld/2)/sqrt(alphasq))),
               'one_peak': 0.5*2/(exp(x-203./15.*40./Ld)+exp(-x+203./15.*40./Ld)),
               'proper_peak': 0.5*2/(exp(x-Ld/4)+exp(-x+Ld/4)),
               'flat': Constant(2*pi**2/(9*40**2)),
               'fast_flat': Constant(0.1),
               'coshes': Constant(2000)*cosh((2000**0.5/2)*(x-0.75))**(-2)+Constant(1000)*cosh(1000**0.5/2*(x-0.25))**(-2),
               'd_peakon':exp(-sqrt((x-Ld/2)**2 + epsilon * deltax ** 2) / sqrt(alphasq)),
               'zero': Constant(0.0),
               'two_peakons': conditional(x < Ld/4, exp((x-Ld/4)/sqrt(alphasq)) - exp(-(x+Ld/4)/sqrt(alphasq)),
                                          conditional(x < 3*Ld/4, exp(-(x-Ld/4)/sqrt(alphasq)) - exp((x-3*Ld/4)/sqrt(alphasq)),
                                                      exp((x-5*Ld/4)/sqrt(alphasq)) - exp(-(x-3*Ld/4)/sqrt(alphasq)))),
               'twin_peakons': conditional(x < Ld/4, exp((x-Ld/4)/sqrt(alphasq)) + 0.5* exp((x-Ld/2)/sqrt(alphasq)),
                                           conditional(x < Ld/2, exp(-(x-Ld/4)/sqrt(alphasq)) + 0.5* exp((x-Ld/2)/sqrt(alphasq)),
                                                       conditional(x < 3*Ld/4, exp(-(x-Ld/4)/sqrt(alphasq)) + 0.5 * exp(-(x-Ld/2)/sqrt(alphasq)),
                                                                   exp((x-5*Ld/4)/sqrt(alphasq)) + 0.5 * exp(-(x-Ld/2)/sqrt(alphasq))))),
               'periodic_peakon': conditional(x < Ld/2, exp((x-Ld/2)/sqrt(alphasq)) + exp(-Ld)*exp(-(x-Ld/2)/sqrt(alphasq)),
                                                        exp(-(x-Ld/2)/sqrt(alphasq)) + exp(-Ld/sqrt(alphasq))*exp((x-Ld/2)/sqrt(alphasq))),
               'cos_bell':conditional(x < Ld/4, (cos(pi*(x-Ld/8)/(2*Ld/8)))**2, 0.0),
               'antisymmetric': 1/(exp((x-Ld/4)/Ld)+exp((-x+Ld/4)/Ld)) - 1/(exp((Ld-x-Ld/4)/Ld)+exp((Ld+x+Ld/4)/Ld))}

    ic_expr = ic_dict[ic]

    if prognostic_variables.scheme in ['upwind', 'LASCH']:

        VCG5 = FunctionSpace(mesh, "CG", 5)
        smooth_condition = Function(VCG5).interpolate(ic_expr)
        prognostic_variables.u.project(as_vector([smooth_condition]))

        # need to find initial m by solving helmholtz problem
        CG1 = FunctionSpace(mesh, "CG", 1)
        u0 = prognostic_variables.u
        p = TestFunction(CG1)
        m_CG = Function(CG1)
        ones = Function(prognostic_variables.Vu).project(as_vector([Constant(1.)]))

        Lm = (p*m_CG - p*dot(ones,u0) - alphasq*p.dx(0)*dot(ones,u0.dx(0)))*dx
        mprob0 = NonlinearVariationalProblem(Lm, m_CG)
        msolver0 = NonlinearVariationalSolver(mprob0, solver_parameters={'ksp_type': 'preonly',
                                                                         'pc_type': 'lu'})
        msolver0.solve()
        prognostic_variables.m.interpolate(m_CG)

        if prognostic_variables.scheme == 'LASCH':
            prognostic_variables.Eu.assign(prognostic_variables.u)
            prognostic_variables.Em.assign(prognostic_variables.m)

    elif prognostic_variables.scheme in ('conforming', 'hydrodynamic', 'test', 'LASCH_hydrodynamic','LASCH_hydrodynamic_m', 'no_gradient'):
        if ic == 'peakon':
            Vu = prognostic_variables.Vu
            # delta = Function(Vu)
            # middle_index = int(len(delta.dat.data[:]) / 2)
            # delta.dat.data[middle_index] = 1
            # u0 = prognostic_variables.u
            # phi = TestFunction(Vu)
            #
            # eqn = phi * u0 * dx + alphasq * phi.dx(0) * u0.dx(0) * dx - phi * delta * dx
            # prob = NonlinearVariationalProblem(eqn, u0)
            # solver = NonlinearVariationalSolver(prob)
            # solver.solve()
            # W = MixedFunctionSpace((Vu, Vu))
            # psi, phi = TestFunctions(W)
            # w = Function(W)
            # u, F = w.split()
            # u.interpolate(ic_expr)
            # u, F = split(w)
            #
            # eqn = (psi * u * dx - psi * (0.5 * u * u + F) * dx
            #        + phi * F * dx + alphasq * phi.dx(0) * F.dx(0) * dx
            #        - phi * u * u * dx - 0.5 * alphasq * phi * u.dx(0) * u.dx(0) * dx)
            #
            # u, F = w.split()
            #
            # prob = NonlinearVariationalProblem(eqn, w)
            # solver = NonlinearVariationalSolver(prob)
            # solver.solve()
            # prognostic_variables.u.assign(u)
            prognostic_variables.u.project(ic_expr)
            # prognostic_variables.u.interpolate(ic_expr)
        else:
            VCG5 = FunctionSpace(mesh, "CG", 5)
            smooth_condition = Function(VCG5).interpolate(ic_expr)
            prognostic_variables.u.project(smooth_condition)

        if prognostic_variables.scheme in ['LASCH_hydrodynamic', 'LASCH_hydrodynamic_m']:
            prognostic_variables.Eu.assign(prognostic_variables.u)

    else:
        raise NotImplementedError('Other schemes not yet implemented.')
