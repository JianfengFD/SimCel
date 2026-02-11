#!/usr/bin/env python3
"""
Visualize MeditRBC simulation output (Corrected Version).

Features:
  - Correctly handles signed edges in E_F (Surface Evolver topology).
  - Calculates correct tangent frames (T1, T2) for curved surfaces.
  - Uses PyVista for high-quality 3D rendering (arrows, scalar fields).

Usage:
  1) Shape only:         python visualize.py shape.dat
  2) Shape + phi:        python visualize.py shape.dat --phi phi.dat
  3) Shape + phi + Q:    python visualize.py shape.dat --phi phi.dat --Q Q.dat
"""

import argparse
import numpy as np
import pyvista as pv

# =============================================================
#  File readers
# =============================================================

def read_shape(filename):
    """
    Read Surface Evolver .dat file.
    Returns:
        rv_dict : dict vertex_index -> np.array([x,y,z])
        V_E     : dict edge_index   -> (v1, v2)  [1-based indices]
        E_F     : dict face_index   -> [e1, e2, e3] (signed edges)
    """
    rv = {}
    V_E = {}
    E_F = {}
    section = None

    with open(filename, 'r') as f:
        for line in f:
            s = line.strip()
            # Skip comments and empty lines
            if not s or s.startswith('//'):
                # Heuristic: reset section if we see a comment block separator
                if 'vertices' in s.lower(): section = 'vertices'
                elif 'edges' in s.lower(): section = 'edges'
                elif 'faces' in s.lower(): section = 'faces'
                continue
            
            # Identify sections by keywords
            if s.lower() == 'vertices':
                section = 'vertices'; continue
            elif s.lower() == 'edges':
                section = 'edges'; continue
            elif s.lower() == 'faces':
                section = 'faces'; continue

            parts = s.split()
            
            try:
                if section == 'vertices' and len(parts) >= 4:
                    idx = int(parts[0])
                    rv[idx] = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
                    
                elif section == 'edges' and len(parts) >= 3:
                    idx = int(parts[0])
                    # Store edge vertices. Note: Surface Evolver edges are 1-based.
                    V_E[idx] = (int(parts[1]), int(parts[2]))
                    
                elif section == 'faces' and len(parts) >= 4:
                    idx = int(parts[0])
                    # Read edges. Note: Surface Evolver faces can represent polygons, 
                    # but usually triangles (3 edges). We read all provided edge indices.
                    # Stop reading if we hit a non-integer (like 'color')
                    edges = []
                    for p in parts[1:]:
                        try:
                            val = int(p)
                            edges.append(val)
                        except ValueError:
                            break
                    E_F[idx] = edges
            except ValueError:
                continue

    return rv, V_E, E_F


def build_mesh_topology(rv_dict, V_E, E_F):
    """
    Convert dicts to numpy arrays respecting topology (signed edges).
    
    CRITICAL FIX: 
    Surface Evolver uses signed edges in faces.
    Edge k: v1 -> v2
    Face uses  k: Traverse v1 -> v2 (Start vertex is v1)
    Face uses -k: Traverse v2 -> v1 (Start vertex is v2)
    
    Returns:
        rv  : (N_V+1, 3) vertex positions
        tri : (N_F, 3)   triangle vertex indices (1-based)
        N_V : max vertex index
    """
    if not rv_dict:
        raise ValueError("No vertices found in shape file.")
        
    N_V = max(rv_dict.keys())
    # Create array with size N_V + 1 to accomodate 1-based indexing directly
    rv = np.zeros((N_V + 1, 3))
    for idx, pos in rv_dict.items():
        rv[idx] = pos

    N_F = len(E_F)
    tri = np.zeros((N_F, 3), dtype=int)

    # Sort keys to ensure deterministic order
    sorted_faces = sorted(E_F.keys())

    for i, fi in enumerate(sorted_faces):
        signed_edges = E_F[fi]
        # We assume the face is a triangle or we take the first 3 vertices
        face_verts = []
        
        for e_signed in signed_edges:
            e_abs = abs(e_signed)
            if e_abs not in V_E:
                continue # Should not happen in valid files
                
            v1, v2 = V_E[e_abs]
            
            # Logic: If edge is positive, we are at v1. If negative, we are at v2.
            # This preserves the winding order v_a -> v_b -> v_c
            if e_signed > 0:
                face_verts.append(v1)
            else:
                face_verts.append(v2)
        
        # Take the first 3 vertices to form a triangle
        if len(face_verts) >= 3:
            tri[i] = face_verts[:3]
        else:
            # Handle degenerate case if necessary
            pass

    return rv, tri, N_V


def read_phi(filename):
    """Read phi file: header line, then vertex_index phi_value"""
    phi = {}
    try:
        with open(filename, 'r') as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith('//'): continue
                parts = s.split()
                if len(parts) >= 2:
                    try:
                        idx = int(parts[0])
                        val = float(parts[1])
                        phi[idx] = val
                    except ValueError:
                        continue
    except FileNotFoundError:
        print(f"Warning: Phi file {filename} not found.")
    return phi


def read_Q(filename):
    """Read Q-tensor file: header, then vertex_index q1 q2 S D"""
    Q = {}
    try:
        with open(filename, 'r') as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith('//'): continue
                parts = s.split()
                if len(parts) >= 4:
                    try:
                        idx = int(parts[0])
                        q1 = float(parts[1])
                        q2 = float(parts[2])
                        S  = float(parts[3])
                        D  = float(parts[4]) if len(parts) >= 5 else 0.0
                        Q[idx] = (q1, q2, S, D)
                    except ValueError:
                        continue
    except FileNotFoundError:
        print(f"Warning: Q file {filename} not found.")
    return Q


# =============================================================
#  Geometry helpers (Corrected)
# =============================================================

def compute_vertex_normals(rv, tri, N_V):
    """
    Compute area-weighted vertex normals.
    Requires 'tri' to have correct winding order (CCW).
    """
    normals = np.zeros((N_V + 1, 3))
    
    # Vectorized calculation is faster
    # tri indices are 1-based, consistent with rv array structure
    v0 = rv[tri[:, 0]]
    v1 = rv[tri[:, 1]]
    v2 = rv[tri[:, 2]]
    
    e1 = v1 - v0
    e2 = v2 - v0
    
    # Face normals (un-normalized, length = 2*area)
    face_normals = np.cross(e1, e2)
    
    # Accumulate onto vertices
    for i in range(3):
        # np.add.at handles repeated indices correctly
        np.add.at(normals, tri[:, i], face_normals)

    # Normalize
    norm = np.linalg.norm(normals, axis=1, keepdims=True)
    # Avoid division by zero
    norm[norm < 1e-15] = 1.0
    normals /= norm
    
    return normals


def compute_tangent_frame(rv, tri, normals, V_E, N_V):
    """
    Build (T1, T2) tangent frame at each vertex.
    T1 = normalized projection of an adjacent edge.
    T2 = n x T1.
    """
    T1 = np.zeros((N_V + 1, 3))
    T2 = np.zeros((N_V + 1, 3))

    # Pre-build an adjacency map: vertex -> one neighbor
    adj = {}
    for (v_start, v_end) in V_E.values():
        if v_start not in adj: adj[v_start] = v_end
        if v_end not in adj:   adj[v_end] = v_start

    for i in range(1, N_V + 1):
        if i not in adj: continue
        
        n = normals[i]
        neighbor = adj[i]
        
        # Edge vector
        edge = rv[neighbor] - rv[i]
        
        # Project onto tangent plane: e_tan = e - (e.n)n
        # This ensures T1 is strictly orthogonal to n
        proj = edge - np.dot(edge, n) * n
        
        length = np.linalg.norm(proj)
        if length < 1e-15:
            # Fallback for degenerate geometry (rare)
            # Pick arbitrary vector
            arb = np.array([1.0, 0.0, 0.0])
            if abs(n[0]) > 0.9: arb = np.array([0.0, 1.0, 0.0])
            proj = np.cross(n, arb)
            length = np.linalg.norm(proj)
        
        t1 = proj / length
        t2 = np.cross(n, t1)
        
        T1[i] = t1
        T2[i] = t2

    return T1, T2


# =============================================================
#  Visualization (PyVista)
# =============================================================

def visualize(rv, tri, N_V, phi_dict=None, Q_dict=None, V_E=None, args=None):
    """
    Main visualization routine using PyVista.
    """
    print("Building PyVista mesh...")
    
    # 1. Prepare Mesh Data
    # PyVista expects faces as [n_verts, v1, v2, v3, n_verts, v1, ...]
    # tri is (N_F, 3). We prepend a column of 3s.
    n_faces = tri.shape[0]
    padding = np.full((n_faces, 1), 3, dtype=int)
    faces_formatted = np.hstack((padding, tri)).flatten()
    
    # Create PolyData object
    mesh = pv.PolyData(rv, faces_formatted)

    # 2. Add Scalar Field (Phi or S)
    scalar_name = None
    cmap = 'coolwarm'
    
    # Initialize full array (including index 0)
    point_data = np.zeros(N_V + 1)
    
    if phi_dict:
        print("Mapping Phi field...")
        # Fill data
        keys = np.array(list(phi_dict.keys()))
        vals = np.array(list(phi_dict.values()))
        # Filter keys within range
        valid = keys <= N_V
        point_data[keys[valid]] = vals[valid]
        
        scalar_name = "Phi"
        cmap = "coolwarm"
        
    elif Q_dict:
        print("Mapping Order Parameter S...")
        for i in range(1, N_V + 1):
            if i in Q_dict:
                point_data[i] = Q_dict[i][2] # S
        scalar_name = "S"
        cmap = "viridis"
    else:
        # Default color for shape-only
        scalar_name = "Elevation"
        point_data = rv[:, 2] # Color by Z
        cmap = "bone"

    mesh.point_data[scalar_name] = point_data

    # 3. Add Director Field (Arrows) if Q is present
    arrows = None
    if Q_dict and V_E:
        print("Calculating director field...")
        # Compute Frames
        normals = compute_vertex_normals(rv, tri, N_V)
        T1, T2 = compute_tangent_frame(rv, tri, normals, V_E, N_V)
        
        arrow_pos = []
        arrow_vec = []
        
        # Subsamping for visibility
        n_arrows_target = getattr(args, 'arrows', 500)
        vert_ids = sorted(list(Q_dict.keys()))
        step = max(1, len(vert_ids) // n_arrows_target)
        
        # Auto-scaling
        bbox = rv[1:].max(axis=0) - rv[1:].min(axis=0)
        scale_len = np.linalg.norm(bbox) * 0.03
        
        for idx in vert_ids[::step]:
            if idx > N_V: continue
            
            q1, q2, S, D = Q_dict[idx]
            
            # Reconstruct director n = cos(alpha)T1 + sin(alpha)T2
            # q1 = S * cos(2*alpha), q2 = S * sin(2*alpha)
            # alpha = 0.5 * atan2(q2, q1)
            alpha = 0.5 * np.arctan2(q2, q1)
            
            u = np.cos(alpha) * T1[idx] + np.sin(alpha) * T2[idx]
            
            # Visualization length proportional to Order Parameter S
            vec_len = scale_len * S
            
            # Offset slightly along normal to avoid Z-fighting
            pos = rv[idx] + normals[idx] * (scale_len * 0.2)
            
            arrow_pos.append(pos)
            arrow_vec.append(u * vec_len)
            
        if arrow_pos:
            # Create a separate PolyData for arrows
            pdata_arrows = pv.PolyData(np.array(arrow_pos))
            pdata_arrows.point_data["vec"] = np.array(arrow_vec)
            # Use glyph filter to create arrows
            arrows = pdata_arrows.glyph(orient="vec", scale=False, geom=pv.Arrow(tip_length=0.25, tip_radius=0.1))

    # 4. Plotting
    pl = pv.Plotter()
    pl.set_background('white')
    
    # Add surface
    pl.add_mesh(mesh, scalars=scalar_name, cmap=cmap, 
                show_edges=False, smooth_shading=True, specular=0.3, label='Surface')
    
    # Add arrows
    if arrows:
        pl.add_mesh(arrows, color='black', label='Director')

    # Setup camera
    pl.add_axes()
    pl.camera_position = 'iso'
    
    print("Displaying...")
    # Save if requested
    if args and args.output and args.no_show:
        pl.screenshot(args.output)
        print(f"Saved to {args.output}")
    else:
        pl.show()


# =============================================================
#  Main
# =============================================================

def main():
    parser = argparse.ArgumentParser(
        description='Visualize MeditRBC simulation (PyVista Version)',
        epilog="Example: python visualize.py shape.dat --phi phi.dat --Q Q.dat")
    
    parser.add_argument('shape', help='Shape file (.dat)')
    parser.add_argument('--phi', help='Concentration file (.dat)')
    parser.add_argument('--Q', help='Q-tensor file (.dat)')
    parser.add_argument('-o', '--output', help='Save screenshot to file')
    parser.add_argument('--no-show', action='store_true', help='Do not open window')
    parser.add_argument('--arrows', type=int, default=500, help='Number of arrows')
    
    args = parser.parse_args()

    # 1. Read Shape
    print(f"Reading shape: {args.shape}")
    rv_dict, V_E, E_F = read_shape(args.shape)
    
    # 2. Build Topology (Corrected)
    rv, tri, N_V = build_mesh_topology(rv_dict, V_E, E_F)
    print(f"  Vertices: {N_V}, Faces: {len(tri)}")

    # 3. Read Fields
    phi_dict = None
    if args.phi:
        print(f"Reading phi: {args.phi}")
        phi_dict = read_phi(args.phi)

    Q_dict = None
    if args.Q:
        print(f"Reading Q: {args.Q}")
        Q_dict = read_Q(args.Q)

    # 4. Visualize
    visualize(rv, tri, N_V, phi_dict, Q_dict, V_E, args)

if __name__ == '__main__':
    main()