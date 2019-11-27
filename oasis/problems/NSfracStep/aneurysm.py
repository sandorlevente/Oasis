__author__ = "Mikael Mortensen <mikaem@math.uio.no>"
__date__ = "2013-06-25"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__ = "GNU Lesser GPL version 3 or any later version"

from ..NSfracStep import *


# Override some problem specific parameters
def problem_parameters(NS_parameters, commandline_kwargs, NS_expressions, **NS_namespace):
    NS_parameters.update(
        folder="aneurysm_results",

        nu=0.0001,
        T=10,
        dt=0.1,

 
        velocity_degree=1,   # polynomial degrees
        pressure_degree=1,    

        # Some discretization options
        # Use Adams Bashforth projection as first estimate for pressure on new timestep
        AB_projection_pressure=False,
        solver="IPCS_ABCN",  # "IPCS_ABCN", "IPCS_ABE", "IPCS", "Chorin", "BDFPC", "BDFPC_Fast"

        # Parameters used to tweak solver
        max_iter=2,                 # Number of inner pressure velocity iterations on timestep
        max_error=1e-6,             # Tolerance for inner iterations (pressure velocity iterations)
        iters_on_first_timestep=2,  # Number of iterations on first timestep
        print_intermediate_info=0, #Summary of timings
        print_velocity_pressure_convergence=True,

        # Parameters used to tweak output
        plot_interval=0,
        checkpoint=10,              # Overwrite solution in Checkpoint folder each checkpoint
        save_step=10,               # Store solution each save_step
        restart_folder=None,        # If restarting solution, set the folder holding the solution to start from here
        output_timeseries_as_vector=True,  # Store velocity as vector in Timeseries


        use_krylov_solvers=True,
        krylov_solvers=dict(
            monitor_convergence=True,
            report=True,
            error_on_nonconvergence=False,
            nonzero_initial_guess=True,
            maximum_iterations=10,
            relative_tolerance=1e-8,
            absolute_tolerance=1e-8)
    )

    NS_expressions.update(dict(p_in=Expression("10*sin(pi*t)", t=0., degree=1)))


mesh = Mesh("mesh/aneurysm.xml")
print(mesh.geometric_dimension(),"D mesh imported")
print("Number of cells: ",mesh.num_cells())
print("Maximum cell size: ",mesh.hmax())
print("Minimum cell size: ",mesh.hmin())
info_green('Mesh imported')

boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1, mesh.domains())
# write out patch IDs
File('boundaries.pvd') << boundaries
info_green('BCs loaded')


def create_bcs(V, Q, sys_comp, p_in, **NS_namespace):

    bcs = dict((ui, []) for ui in sys_comp)
    wall = DirichletBC(V, 0.0, boundaries, 0)
    p_inlet = DirichletBC(Q, 10, boundaries, 2)
    bcp = DirichletBC(Q, 0.0, boundaries, 1)
    bcp2 = DirichletBC(Q, 0.0, boundaries, 3)
    # Boundary conditions for each u,v,w velocity component
    bcs['u0']=[wall]
    bcs['u1']=[wall]
    bcs['u2']=[wall]
    # BC for pressure
    bcs['p']=[bcp,bcp2, p_inlet]
    return bcs

##ez is mi a tököm??

"""def pre_solve_hook(mesh, velocity_degree, u_,
                   AssignedVectorFunction, **NS_namespace):
    return dict(uv=AssignedVectorFunction(u_))
"""

def start_timestep_hook(t, p_in, **NS_namespace):
    p_in.t = t



# ez mi a szarra kell?
def initialize(x_1, x_2,bcs, **NS_namespace):
    for ui in x_2:
        [bc.apply(x_1[ui]) for bc in bcs[ui]]
        [bc.apply(x_2[ui]) for bc in bcs[ui]]


def temporal_hook(tstep,t, u_, plot_interval, **NS_namespace):
    print("Time= ",t, " at simulation timestep =", tstep)
    #if tstep % plot_interval == 0:
        
        # print(type(u_))
        #print(u_components)
    #    print('u max:', u_[0].max()) 
