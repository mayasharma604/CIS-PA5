import numpy as np

'''
Created on October 4, 2025
Author: Maya Sharma
FrameTransform: implements the registration algorithm in 9/9 lecture slides 17
params: source and target point clouds (N x 3) - N points with xyz coordinates
returns: R and t: rotation matrix and translation vector that transforms points A to B
summary: centers both point clouds, computes covariance matrix, performs SVD on H, and constructs R. 
makes sure that R is a rotation matrix, not reflection. Lastly, computes translation.
'''

def register_points(A, B):
    
    assert A.shape == B.shape 
    # make sure both are 3D
    
    # get centroids
    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)
    
    # center points
    AA = A - centroid_A
    BB = B - centroid_B
    
    # covariance of matrix H
    H = AA.T @ BB
    
    # calculate SVD of H
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    
    # check det(R) = 1 is true, otherwise correct for relfection
    if np.linalg.det(R) < 0:
        Vt[2, :] *= -1
        R = Vt.T @ U.T
    
    # translation
    t = centroid_B - R @ centroid_A
    
    return R, t

'''
Created on October 2, 2025
Author: Maya Sharma
params: rotation and translation matris as well as coord points
Returns: transformed points (N x 3)
Summary: apply basic transformation to frame
'''
def apply_transform(R, t, pts):
    return (R @ pts.T).T + t


'''
Updated on December 8, 2025
Author: Maya Sharma
params: rotation and translation matris 
Returns: inversely transformed points 
Summary: compute inverse of transformation F = (R,t)
'''
def apply_inverse_transform(R, t):
    R_inv = R.T 
    t_inv = -t @ R_inv 
    
    
    return R_inv, t_inv

'''
Created on November 9, 2025
Author: Maya Sharma
params: file name 
Returns: markers (nx3) and the tip (1x3)
Summary: reads body definition file - handles both space and comma separated formats
'''
def read_body(filename):
    
    with open(filename, 'r') as f:
        lines = f.readlines()

    header_line = lines[0].strip()
    
    if ',' in header_line:
        header = [x.strip() for x in header_line.split(',')]
    else:
        header = [x for x in header_line.split() if x]
    
    Nmarkers = int(header[0])
    
    markers = []
    for i in range(1, Nmarkers + 1):
        line = lines[i].strip()
        if ',' in line:
            coords = [x.strip() for x in line.split(',')]
        else:
            coords = [x for x in line.split() if x]
        marker_coords = [float(x) for x in coords]
        markers.append(marker_coords)
    markers = np.array(markers)
    
    tip_line = lines[Nmarkers + 1].strip()
    if ',' in tip_line:
        tip_coords = [x.strip() for x in tip_line.split(',')]
    else:
        tip_coords = [x for x in tip_line.split() if x]
    tip_coords = [float(x) for x in tip_coords]
    tip = np.array(tip_coords)
    
    return markers, tip

'''
Created on November 6, 2025
Author: Maya Sharma
params: filename
Returns: vertices, triangles, neighbours
Summary: reads mesh file
'''
def read_mesh(filename):
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    n_vertices = int(lines[0].strip())
    
    vertices = []
    for i in range(1, n_vertices + 1):
        coords = [x for x in lines[i].strip().split() if x]
        v_coords = [float(x) for x in coords]
        vertices.append(v_coords)
    vertices = np.array(vertices)
    
    n_triangles = int(lines[n_vertices + 1].strip())
    
    triangles = []
    neighbours = []
    for i in range(n_vertices + 2, n_vertices + 2 + n_triangles):
        indices = [x for x in lines[i].strip().split() if x]
        indices = [int(x) for x in indices]
        triangles.append(indices[:3])
        neighbours.append(indices[3:])
    triangles = np.array(triangles)
    neighbours = np.array(neighbours)
    
    return vertices, triangles, neighbours

'''
Updated on December 7, 2025
Author: Maya Sharma
params: filename
Returns: frames (which is a list of nx3 arrays), Ns, Nsamps
Summary: reads sample readings file
'''
def read_sample_readings(filename):
    """Read sample readings file - returns 4 values."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Parse header
    header_line = lines[0].strip()
    header_parts = [x.strip() for x in header_line.split(',')]
    
    Ns = int(header_parts[0])
    Nsamps = int(header_parts[1])
    
    # Extract N_modes from last part: "filename N_modes"
    filename_and_modes = header_parts[2]
    N_modes = int(filename_and_modes.split()[-1])
    
    frames = []
    idx = 1
    for frame in range(Nsamps):
        records = []
        for i in range(Ns):
            line = lines[idx].strip()
            coords = [float(x.strip()) for x in line.split(',') if x.strip()]
            records.append(coords)
            idx += 1
        frames.append(np.array(records))
    
    # RETURN 4 VALUES, not 3
    return frames, Ns, Nsamps, N_modes
