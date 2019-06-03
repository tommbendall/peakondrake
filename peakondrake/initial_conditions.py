from firedrake import (dx, conditional, exp, as_vector, dot, pi,
                       Function, NonlinearVariationalSolver,
                       TestFunction, NonlinearVariationalProblem,
                       SpatialCoordinate, Constant, FunctionSpace, cosh)

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

    ic_dict = {'two_peaks': (0.2*2/(exp(x-403./15.) + exp(-x+403./15.))
                             + 0.5*2/(exp(x-203./15.)+exp(-x+203./15.))),
               'gaussian': 0.5*exp(-((x-10.)/2.)**2),
               'gaussian_narrow': 0.5*exp(-((x-10.)/1.)**2),
               'gaussian_wide': 0.5*exp(-((x-10.)/3.)**2),
               'peakon': conditional(x < 20., exp((x-20.)/1.), exp(-(x-20.)/1.)),
               'one_peak': 0.5*2/(exp(x-203./15.)+exp(-x+203./15.)),
               'flat': Constant(2*pi**2/(9*40**2)),
               'coshes': Constant(2000)*cosh((2000**0.5/2)*(x-0.75))**(-2)+Constant(1000)*cosh(1000**0.5/2*(x-0.25))**(-2)}

    ic_expr = ic_dict[ic]

    if prognostic_variables.scheme == 'upwind':

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

    elif prognostic_variables.scheme == 'conforming':
        VCG5 = FunctionSpace(mesh, "CG", 5)
        smooth_condition = Function(VCG5).interpolate(ic_expr)
        prognostic_variables.u.project(smooth_condition)

    else:
        raise NotImplementedError('Other schemes not yet implemented.')
