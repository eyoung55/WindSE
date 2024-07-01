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
name = "velocity_probe"

### Run this code to see which names are already taken ###
# import windse.objective_functions as obj_funcs
# print(obj_funcs.objective_functions.keys())
# print(obj_funcs.objective_kwargs)
# exit()

### Set default keyword argument values ###
# These must be a dictionary and will be passed in via the kwargs.
# Leave empty if no argument are needed. 
keyword_defaults = {}

### Define objective function
def objective(solver, inflow_angle = 0.0, first_call=False, **kwargs):
    '''
    This uses a dolfin function eval to measure the 
    velocity at a point located at the rotor hub, i.e., the 
    (x, y, hub_height) location, where (x, y) is the location of 
    turbine and hub_height is the hub height (modified to reflect 
    distance from ground if terrain is present).
    '''

    velocity_array = np.zeros((solver.problem.farm.numturbs, 9))

    print('Calculating Velocity at Each Turbine')

    for k in range(solver.problem.farm.numturbs):
        # Get the (x, y, z) position of this turbine centroid
        mx =  solver.problem.farm.turbines[k].mx
        my =  solver.problem.farm.turbines[k].my
        mz =  solver.problem.farm.turbines[k].mz

        # Use the (x, y, z) position to get the velocity
        velocity_at_point = solver.problem.u_k(mx, my, mz)

        u_x = velocity_at_point[0]
        u_y = velocity_at_point[1]
        u_z = velocity_at_point[2]
        u_mag = np.sqrt(u_x * u_x + u_y * u_y + u_z * u_z)

        # Get the yaw for writing the file
        yaw = solver.problem.farm.turbines[k].myaw

        # Store the turbine id, x, y, z, yaw, and velocity 
        # on the kth row for turbine #k
        velocity_array[k, :] = [k, mx, my, mz, yaw, u_x, u_y, u_z, u_mag]

    folder_string = solver.params.folder+"data/"
    np.savetxt('%sturbine_velocities.csv' % (folder_string),
        velocity_array,
        fmt='%.6e',
        header='Turbine ID (#), X-Location (m), Y-Location (m), Z-Location (m), Yaw (Rad), X-Velocity (m/s), Y-Velocity (m/s), Z-Velocity (m/s), Velocity Magnitude (m/s)',
        delimiter=',')

    # This will be the average of the 8th column, that is,
    # the average of all the point-wise velocity magnitudes calculated above
    J = np.mean(velocity_array[:, 8])

    print('Objective Value: ', float(J))

    return J
