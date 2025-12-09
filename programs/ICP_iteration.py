import numpy as np
from utility_functions import *
from ICP_algo import *

'''
Created on December 2, 2025
Author: Maya Sharma
Parameters: Body A/B markers and tip, sample readings data
Returns: Array of d_k points (N_samps x 3)
Summary: Pre-calculates the pointer tip position d_k in the Body B coordinate system for all sample frames.
This value is constant throughout the ICP iterations.
'''
def pre_calculate_dks(A_markers, A_tip, B_markers, B_tip, frames, Nsamps, N_A, N_B):
    d_k_points = []
    

    for k in range(Nsamps):
        frame_data = frames[k]
        
        a_markers_tracker = frame_data[:N_A]
        b_markers_tracker = frame_data[N_A:N_A+N_B]
        
        # calculate rigid bodies
        R_A, t_A = register_points(A_markers, a_markers_tracker)
        R_B, t_B = register_points(B_markers, b_markers_tracker)
        
        R_B_inv, t_B_inv = apply_inverse_transform(R_B, t_B)
        
        # Transform A_tip from Body A to Tracker frame
        A_tip_transformed = apply_transform(R_A, t_A, A_tip.reshape(1, -1))[0]
        # Transform result from Tracker frame to Body B frame (d_k)
        d_k = apply_transform(R_B_inv, t_B_inv, A_tip_transformed.reshape(1, -1))[0]
        
        d_k_points.append(d_k)
        
    return np.array(d_k_points)

'''
Created on December 2, 2025
Author: Maya Sharma
Parameters: takes in bodyA_file, bodyB_file, the mesh file, sample readings file, and the output file,
            plus max_iterations and a tolerance for convergence.
Returns: The final s_k and c_k points
Summary: main function for the complete ICP algorithm (PA#4)
'''

def solve_pa4(bodyA_file, bodyB_file, mesh_file, sample_readings_file,
              output_file, max_iterations=50, tolerance=1e-5):

    # load mesh and rigid body definitions
    vertices, triangles, _ = read_mesh(mesh_file)
    A_markers, A_tip = read_body(bodyA_file)
    B_markers, B_tip = read_body(bodyB_file)
    frames, _, Nsamps = read_sample_readings(sample_readings_file)

    N_A = len(A_markers)
    N_B = len(B_markers)

    # precalculate all d_k 
    d_k_points = pre_calculate_dks(
        A_markers, A_tip,
        B_markers, B_tip,
        frames, Nsamps, N_A, N_B
    )

    # Initialize F_reg = identity
    R_reg = np.eye(3)
    t_reg = np.zeros(3)

    prev_mean_dist = float('inf')

    print("Starting ICP iterations...")

    # These will store the last iteration s_k and c_k
    last_s_k_points = None
    last_c_k_points = None

    
    for iteration in range(max_iterations):

        s_k_points = []
        c_k_points = []
        distances = []

        for d_k in d_k_points:
            # s_k = F_reg(d_k)
            s_k = apply_transform(R_reg, t_reg, d_k.reshape(1, -1))[0]

            # Find closest mesh point
            c_k, dist, _ = closest_point_on_mesh(s_k, vertices, triangles)

            s_k_points.append(s_k)
            c_k_points.append(c_k)
            distances.append(dist)

        s_k_points = np.array(s_k_points)
        c_k_points = np.array(c_k_points)

        # save these 
        last_s_k_points = s_k_points.copy()
        last_c_k_points = c_k_points.copy()

        R_new, t_new = register_points(d_k_points, c_k_points)

        mean_dist = np.mean(distances)
        print(f"Iteration {iteration+1}: mean distance = {mean_dist:.6f}")

        # Check convergence
        if abs(prev_mean_dist - mean_dist) < tolerance:
            print(f"Converged at iteration {iteration+1}")
            break

        # Update transform
        R_reg = R_new
        t_reg = t_new
        prev_mean_dist = mean_dist

    else:
        print(f"Reached max iterations ({max_iterations}) without full convergence.")

    final_s = last_s_k_points
    final_c = last_c_k_points

    with open(output_file, "w") as f:
        f.write(f"{Nsamps} {output_file}\n")

        for s_k, c_k in zip(final_s, final_c):
            diff = np.linalg.norm(s_k - c_k)

            f.write(f"{s_k[0]:8.2f} {s_k[1]:8.2f} {s_k[2]:8.2f}    ")
            f.write(f"{c_k[0]:8.2f} {c_k[1]:8.2f} {c_k[2]:8.2f}    ")
            f.write(f"{diff:8.3f}\n")

    return final_s, final_c
