################## HEADER DO NOT EDIT ##################
import os
import __main__

### Get the name of program importing this package ###
if hasattr(__main__,"__file__"):
    main_file = os.path.basename(__main__.__file__)
else:
    main_file = "ipython"
    
### This checks if we are just doing documentation ###
if not main_file in ["sphinx-build", "__main__.py"]:
    from dolfin import *
    from dolfin_adjoint import *
########################################################

### Additional import statements ###
import numpy as np

### Declare Unique name
name = "point_blockage"


### Set default keyword argument values ###
keyword_defaults = {
    "location": (0,0,0),
    "offset_by_mean": False
}


### Define objective function
def objective(solver, inflow_angle = 0.0, first_call=False, **kwargs):
    '''
    This is a simple blockage metric that evaluates the velocity deficit at 
    a single location in the farm.

    Keyword arguments:
        location: where the deficit is evaluated
    '''

    x0 = np.array(kwargs.pop("location"))
    # u_ref, v_ref, w_ref = solver.problem.bd.bc_velocity.split()
    # u,     v,     w     = solver.problem.u_k.split()
    # u = solver.problem.p_k
    # ud = (u(x0)-u_ref(x0))/u_ref(x0)
    # ud = u(x0)
    
    # This doesn't work in parallel (point may not be inside process X's domain)
    # J = solver.problem.up_k(x0)[0]

    # try:
    #     vel_at_point = np.array([solver.problem.up_k(x0)[0]], dtype=np.float64)
    # except:
    #     vel_at_point = np.array([np.nan], dtype=np.float64)

    # gathered_vel_at_point = np.zeros(solver.params.num_procs, dtype=np.float64)
    # solver.params.comm.Allgather(vel_at_point, gathered_vel_at_point)

    x = SpatialCoordinate(solver.problem.dom.mesh)

    delta_x = x[0] - x0[0]
    delta_y = x[1] - x0[1]
    delta_z = x[2] - x0[2]

    distance = (delta_x**2 + delta_y**2 + delta_z**2)/solver.problem.dom.mesh.hmax()
    spherical_gaussian = exp(-pow(distance, 6.0))
    volume = assemble(spherical_gaussian*dx)
    J = assemble(solver.problem.up_k[0]*spherical_gaussian/volume*dx)

    # try:
    #     vel_at_point = solver.problem.up_k(x0)[0]
    # except:
    #     vel_at_point = np.nan

    # gathered_vel_at_point = solver.params.comm.allgather(vel_at_point)

    # Find the first non-NaN value in a vector
    # def get_first_non_nan(data):
    #     for val in data:
    #         if not np.isnan(val):
    #             return val

    #     return val

    # J = get_first_non_nan(gathered_vel_at_point)
    # print('Rank %d holds type ' % (solver.params.rank), type(J))

    # if np.isnan(J):
    #     raise ValueError("Couldn't find point", x0, "anywhere inside the domain.")

    return J
