# deform_registration_correct.py
import numpy as np
import os
from utility_functions import *
from ICP_algo import *
from ICP_iteration import *

def compute_barycentric(point, triangle_vertices):
    A, B, C = triangle_vertices
    v0 = B - A
    v1 = C - A
    v2 = point - A
    
    d00 = np.dot(v0, v0)
    d01 = np.dot(v0, v1)
    d11 = np.dot(v1, v1)
    d20 = np.dot(v2, v0)
    d21 = np.dot(v2, v1)
    
    denom = d00 * d11 - d01 * d01
    if abs(denom) < 1e-12:
        return 1/3, 1/3, 1/3
    
    v = (d11 * d20 - d01 * d21) / denom
    w = (d00 * d21 - d01 * d20) / denom
    u = 1 - v - w
    
    return u, v, w

def read_modes_fixed(filename, max_modes):
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
    
    import re
    n_vertices_match = re.search(r'Nvertices=(\d+)', lines[0])
    n_modes_match = re.search(r'Nmodes=(\d+)', lines[0])
    
    n_vertices = int(n_vertices_match.group(1))
    n_modes_total = int(n_modes_match.group(1))
    n_modes = min(max_modes, n_modes_total)
    
    mean_vertices = []
    modes = [[] for _ in range(n_modes)]
    
    line_idx = 1
    if line_idx < len(lines) and "Mode 0" in lines[line_idx]:
        line_idx += 1
    
    for i in range(n_vertices):
        line = lines[line_idx]
        if ',' in line:
            coords = [float(x.strip()) for x in line.split(',')]
        else:
            coords = [float(x) for x in line.split()]
        mean_vertices.append(coords)
        line_idx += 1
    
    for m in range(n_modes):
        if line_idx < len(lines) and f"Mode {m+1}" in lines[line_idx]:
            line_idx += 1
        
        for i in range(n_vertices):
            line = lines[line_idx]
            if ',' in line:
                coords = [float(x.strip()) for x in line.split(',')]
            else:
                coords = [float(x) for x in line.split()]
            modes[m].append(coords)
            line_idx += 1
    
    return np.array(mean_vertices), [np.array(mode) for mode in modes]

def solve_pa5_correct(bodyA_file, bodyB_file, mesh_file, modes_file, 
                      sample_readings_file, output_file, max_iters=50):
    
    
    
    vertices, triangles, _ = read_mesh(mesh_file)
    A_markers, A_tip = read_body(bodyA_file)
    B_markers, B_tip = read_body(bodyB_file)
    
    frames, Ns, Nsamps, N_modes = read_sample_readings(sample_readings_file)
    
    N_A = len(A_markers)
    N_B = len(B_markers)
    
    
    mean_vertices, mode_vectors = read_modes_fixed(modes_file, N_modes)
    
    
    d_k_points = []
    
    for k in range(Nsamps):
        frame = frames[k]
        a_markers_tracker = frame[:N_A]
        b_markers_tracker = frame[N_A:N_A+N_B]
        
        R_A, t_A = register_points(A_markers, a_markers_tracker)
        R_B, t_B = register_points(B_markers, b_markers_tracker)
        
        R_B_inv, t_B_inv = apply_inverse_transform(R_B, t_B)
        
        A_tip_tracker = apply_transform(R_A, t_A, A_tip.reshape(1, -1))[0]
        
        d_k = apply_transform(R_B_inv, t_B_inv, A_tip_tracker.reshape(1, -1))[0]
        d_k_points.append(d_k)
    
    d_k_points = np.array(d_k_points) 
    
    R_reg = np.eye(3)
    t_reg = np.zeros(3)
    lambdas = np.zeros(N_modes)
    
    
    for iteration in range(max_iters):
        
       
        deformed_vertices = mean_vertices.copy()
        for m in range(N_modes):
            deformed_vertices += lambdas[m] * mode_vectors[m]
        
        R_reg_inv, t_reg_inv = apply_inverse_transform(R_reg, t_reg) 

        s_k_points = apply_transform(R_reg, t_reg, d_k_points)
        
        d_prime_k_points = apply_transform(R_reg_inv, t_reg_inv, s_k_points) 
        
        bary_coords = []
        triangle_indices = []
        
        for d_prime_k in d_prime_k_points:
            m_k_on_deformed, _, tri_idx = closest_point_on_mesh(d_prime_k, deformed_vertices, triangles)
            
            v0_idx, v1_idx, v2_idx = triangles[tri_idx]
            tri_verts_deformed = deformed_vertices[[v0_idx, v1_idx, v2_idx]] 
            
            zeta, xi, psi = compute_barycentric(m_k_on_deformed, tri_verts_deformed)
            
            bary_coords.append((zeta, xi, psi))
            triangle_indices.append(tri_idx)

        
        A = np.zeros((3 * Nsamps, N_modes))
        b = np.zeros(3 * Nsamps) 
        
        for k in range(Nsamps):
            d_prime_k = d_prime_k_points[k]
            tri_idx = triangle_indices[k]
            zeta, xi, psi = bary_coords[k]
            v0_idx, v1_idx, v2_idx = triangles[tri_idx]

            m0_v0, m0_v1, m0_v2 = mean_vertices[[v0_idx, v1_idx, v2_idx]]
            q0_k = zeta * m0_v0 + xi * m0_v1 + psi * m0_v2
            
            b[3*k:3*k+3] = d_prime_k - q0_k
            
            for m in range(N_modes):
                mv0, mv1, mv2 = mode_vectors[m][[v0_idx, v1_idx, v2_idx]]
                a_mk = zeta * mv0 + xi * mv1 + psi * mv2 # q_m,k
                A[3*k:3*k+3, m] = a_mk
                
        lambda_new, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
        lambdas = lambda_new
        
        deformed_vertices_new = mean_vertices.copy()
        for m in range(N_modes):
            deformed_vertices_new += lambdas[m] * mode_vectors[m]
            
        m_k_for_Freg = [] 
        for d_prime_k in d_prime_k_points:
            closest_pt_final, _, _ = closest_point_on_mesh(d_prime_k, deformed_vertices_new, triangles)
            m_k_for_Freg.append(closest_pt_final)
        m_k_for_Freg = np.array(m_k_for_Freg) 
        
       
        R_new, t_new = register_points(d_k_points, m_k_for_Freg)
        
        s_k_new_final = apply_transform(R_new, t_new, d_k_points)
        c_k_new_final = apply_transform(R_new, t_new, m_k_for_Freg)
        
        mean_error = np.mean(np.linalg.norm(s_k_new_final - c_k_new_final, axis=1))

        R_reg = R_new
        t_reg = t_new
        
        print(f"Iter {iteration+1}: error = {mean_error:.6f}, lambdas = {[f'{l:.2f}' for l in lambdas]}")
        
        if mean_error < 0.001:
            print(f"Converged at iteration {iteration+1}")
            break
    

    final_deformed = mean_vertices.copy()
    for m in range(N_modes):
        final_deformed += lambdas[m] * mode_vectors[m]
    
    s_k_final = apply_transform(R_reg, t_reg, d_k_points)
    
    R_reg_inv_final, t_reg_inv_final = apply_inverse_transform(R_reg, t_reg)
    s_k_final_B = apply_transform(R_reg_inv_final, t_reg_inv_final, s_k_final)
    
    m_k_final_B = []
    for s_k_B in s_k_final_B:
        closest_pt_B, _, _ = closest_point_on_mesh(s_k_B, final_deformed, triangles)
        m_k_final_B.append(closest_pt_B)
    m_k_final_B = np.array(m_k_final_B)
    
    c_k_final = apply_transform(R_reg, t_reg, m_k_final_B)

    final_error = np.mean(np.linalg.norm(s_k_final - c_k_final, axis=1))
    
  
    
    output_name = "PA5-A-Debug-Answer.txt"
    
    with open(output_file, 'w') as f:
        f.write(f"{Nsamps} {output_name} {N_modes}\n")
        
        for m in range(N_modes):
            f.write(f"{lambdas[m]:12.4f}")
        f.write("\n")
        
        for k in range(Nsamps):
            s_k = s_k_final[k]
            c_k = c_k_final[k]
            err = np.linalg.norm(s_k - c_k)
            
            f.write(f"{s_k[0]:8.2f}{s_k[1]:9.2f}{s_k[2]:9.2f}   ")
            f.write(f"{c_k[0]:8.2f}{c_k[1]:9.2f}{c_k[2]:9.2f}  ")
            f.write(f"{err:8.3f}\n")
    
    
    return s_k_final, c_k_final, lambdas