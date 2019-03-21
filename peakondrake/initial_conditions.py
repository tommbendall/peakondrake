from firedrake import (dx, conditional, exp, as_vector, dot,
                       Function, NonlinearVariationalSolver,
                       TestFunction, NonlinearVariationalProblem,
                       SpatialCoordinate, Constant)

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

    if prognostic_variables.scheme == 'upwind':

        if ic == 'two_peaks':
            ic_expr = (0.2*2/(exp(x-403./15.) + exp(-x+403./15.))
                        + 0.5*2/(exp(x-203./15.)+exp(-x+203./15.)))
        elif ic == 'gaussian':
            ic_expr = (0.5*exp(-((x-10.)/2.)**2))
        elif ic == 'gaussian_narrow':
            ic_expr = (0.5*exp(-((x-10.)/1.)**2))
        elif ic == 'gaussian_wide':
            ic_expr = (0.5*exp(-((x-10.)/3.)**2))
        elif ic == 'peakon':
            ic_expr = conditional(x < 20., exp((x-20.)/1.), exp(-(x-20.)/1.))
        elif ic == 'one_peak':
            ic_expr = (0.5*2/(exp(x-203./15.)+exp(-x+203./15.)))
        else:
            raise ValueError('Initial condition not recognised.')

        prognostic_variables.u.project(as_vector([ic_expr]))

        # need to find initial m by solving helmholtz problem
        m0 = Function(prognostic_variables.Vm)
        u0 = prognostic_variables.u
        p = TestFunction(prognostic_variables.Vm)
        ones = Function(prognostic_variables.Vu).project(as_vector([Constant(1.)]))

        Lm = (p*m0 - p*dot(ones,u0) - alphasq*p.dx(0)*dot(ones,u0.dx(0)))*dx
        mprob0 = NonlinearVariationalProblem(Lm, m0)
        msolver0 = NonlinearVariationalSolver(mprob0, solver_parameters={'ksp_type': 'preonly',
                                                                         'pc_type': 'lu'})
        msolver0.solve()
        prognostic_variables.m.assign(m0)

    else:
        raise NotImplementedError('Other schemes not yet implemented.')
