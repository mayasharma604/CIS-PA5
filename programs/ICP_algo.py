import numpy as np
from utility_functions import *


'''
Created on November 6, 2025
Author: Anishka Bhartiya
Parameters: point and triangle_vertices (dimensions 3,3)
Returns: the closest point on triangle to given point
Summary: computes closest point that is on the triangle to given point
'''
def closest_point_on_triangle(point, triangle_vertices):
    # triangle verts
    A, B, C = triangle_vertices
    
    AB = B - A
    AC = C - A
    AP = point - A
    
    d_1 = np.dot(AB, AP)
    d_2 = np.dot(AC, AP)
    
    if d_1 <= 0.0 and d_2 <= 0.0:
        return A
    
    BP = point - B
    
    d3 = np.dot(AB, BP)
    d4 = np.dot(AC, BP)
    if d3 >= 0.0 and d4 <= d3:
        return B
    
    
    CP = point - C
    
    d5 = np.dot(AB, CP)
    d6 = np.dot(AC, CP)
    if d6 >= 0.0 and d5 <= d6:
        return C
    
    vc = d_1 * d4 - d3 * d_2
    
    if vc <= 0.0 and d_1 >= 0.0 and d3 <= 0.0:
        # fraction along AB is
        v = d_1 / (d_1 - d3)
        return A + v * AB
    
    vb = d5 * d_2 - d_1 * d6
    
    if vb <= 0.0 and d_2 >= 0.0 and d6 <= 0.0:
        
        w = d_2 / (d_2 - d6)
        return A + w * AC
    
    va = d3 * d6 - d5 * d4
    if va <= 0.0 and (d4 - d3) >= 0.0 and (d5 - d6) >= 0.0:
        
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        return B + w * (C - B)
    
    denom = 1.0 / (va + vb + vc)
    
    v = vb * denom
    w = vc * denom
    
    return A + v * AB + w * AC

'''
Created November 6, 2025
Author: Anishka Bhartiya
Parameters: takes in the point, vertcies on mesh and triangles on mesh 
Returns: returns the closest point that's on mesh, what the distnce to that point is, and that triangles' index
Summary: this functions finds what the closest point on the mesh to a given point
'''
def closest_point_on_mesh(point, vertices, triangles):
    
    min_distance = float('inf')
    close_pt = None
    closest_triangle_idx = -1
    
    for i, triangle in enumerate(triangles):
        
        # indices in the triangle array are 0-based for numpy indexing
        triangle_vertices = vertices[triangle]
        
        candidate_point = closest_point_on_triangle(point, triangle_vertices)
        distance = np.linalg.norm(point - candidate_point)
        
        if distance < min_distance:
            min_distance = distance
            close_pt = candidate_point
            closest_triangle_idx = i
    
    return close_pt, min_distance, closest_triangle_idx


