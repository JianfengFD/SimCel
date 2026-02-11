#!/usr/bin/env python3
"""
Visualize MeditRBC simulation output.

Three modes:
  1) Shape only:         python visualize.py shape.dat
  2) Shape + phi:        python visualize.py shape.dat --phi phi.dat
  3) Shape + phi + Q:    python visualize.py shape.dat --phi phi.dat --Q Q.dat

Mode 1: Draws the 3D triangulated surface.
Mode 2: Colors each triangle by the concentration field phi(r).
Mode 3: Colors by phi(r) and draws director arrows u_max on the surface.
"""

import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.cm as cm
import matplotlib.colors as mcolors


# =============================================================
#  File readers
# =============================================================

def read_shape(filename):
    """Read Surface Evolver .dat file.
    Returns:
        rv   : dict  vertex_index -> np.array([x,y,z])
        V_E  : dict  edge_index   -> (v1, v2)
        E_F  : dict  face_index   -> (e1, e2, e3)   (signed)
    """
    rv = {}
    V_E = {}
    E_F = {}
    section = None

    with open(filename, 'r') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('//'):
                section = None
                continue
            if s == 'vertices':
                section = 'vertices'; continue
            elif s == 'edges':
                section = 'edges'; continue
            elif s == 'faces':
                section = 'faces'; continue

            parts = s.split()
            if section == 'vertices' and len(parts) >= 4:
                idx = int(parts[0])
                rv[idx] = np.array([float(parts[1]),
                                    float(parts[2]),
                                    float(parts[3])])
            elif section == 'edges' and len(parts) >= 3:
                idx = int(parts[0])
                V_E[idx] = (int(parts[1]), int(parts[2]))
            elif section == 'faces' and len(parts) >= 4:
                idx = int(parts[0])
                E_F[idx] = (int(parts[1]), int(parts[2]), int(parts[3]))

    return rv, V_E, E_F


def build_mesh(rv_dict, V_E, E_F):
    """Convert dicts to numpy arrays.
    Returns:
        rv  : (N_V+1, 3)  vertex positions (1-indexed, rv[0] unused)
        tri : (N_F, 3)    triangle vertex indices (1-indexed)
    """
    N_V = max(rv_dict.keys())
    rv = np.zeros((N_V + 1, 3))
    for idx, pos in rv_dict.items():
        rv[idx] = pos

    N_F = len(E_F)
    tri = np.zeros((N_F, 3), dtype=int)

    for i, fi in enumerate(sorted(E_F.keys())):
        e1, e2, e3 = E_F[fi]
        # Collect vertices from the 3 edges
        vset = set()
        for e in (e1, e2, e3):
            v1, v2 = V_E[abs(e)]
            vset.add(v1)
            vset.add(v2)
        vlist = sorted(vset)
        tri[i, :len(vlist[:3])] = vlist[:3]

    return rv, tri, N_V


def read_phi(filename):
    """Read phi file: header line, then  vertex_index  phi_value"""
    phi = {}
    with open(filename, 'r') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('//'):
                continue
            parts = s.split()
            if len(parts) >= 2:
                try:
                    idx = int(parts[0])
                    val = float(parts[1])
                    phi[idx] = val
                except ValueError:
                    continue
    return phi


def read_Q(filename):
    """Read Q-tensor file: header, then  vertex_index  q1  q2  S  D"""
    Q = {}
    with open(filename, 'r') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('//'):
                continue
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
    return Q


# =============================================================
#  Geometry helpers
# =============================================================

def compute_vertex_normals(rv, tri, N_V):
    """Area-weighted vertex normals."""
    normals = np.zeros((N_V + 1, 3))
    for t in tri:
        v0, v1, v2 = t
        e1 = rv[v1] - rv[v0]
        e2 = rv[v2] - rv[v0]
        fn = np.cross(e1, e2)          # 2 * area * normal
        normals[v0] += fn
        normals[v1] += fn
        normals[v2] += fn

    nrm = np.linalg.norm(normals, axis=1, keepdims=True)
    nrm[nrm < 1e-15] = 1.0
    normals /= nrm
    return normals


def compute_tangent_frame(rv, tri, normals, V_E, N_V):
    """Build (T1, T2) tangent frame at each vertex.

    T1 = projection of an adjacent edge onto the tangent plane, normalised.
    T2 = n x T1.
    """
    T1 = np.zeros((N_V + 1, 3))
    T2 = np.zeros((N_V + 1, 3))

    # Build one-neighbor map from edge list
    neighbor = {}
    for (va, vb) in V_E.values():
        if va not in neighbor:
            neighbor[va] = vb
        if vb not in neighbor:
            neighbor[vb] = va

    for i in range(1, N_V + 1):
        n = normals[i]
        if i not in neighbor:
            continue
        edge = rv[neighbor[i]] - rv[i]
        t1 = edge - np.dot(edge, n) * n
        length = np.linalg.norm(t1)
        if length < 1e-15:
            continue
        t1 /= length
        t2 = np.cross(n, t1)
        T1[i] = t1
        T2[i] = t2

    return T1, T2


def face_scalar(tri, data_dict, default=0.0):
    """Average a per-vertex scalar dict onto faces."""
    vals = np.zeros(len(tri))
    for i, t in enumerate(tri):
        s = 0.0
        for v in t:
            s += data_dict.get(v, default)
        vals[i] = s / 3.0
    return vals


# =============================================================
#  Plotting
# =============================================================

def set_equal_axes(ax, rv, N_V):
    pts = rv[1:N_V + 1]
    lo = pts.min(axis=0)
    hi = pts.max(axis=0)
    mid = 0.5 * (lo + hi)
    half = 0.5 * (hi - lo).max()
    ax.set_xlim(mid[0] - half, mid[0] + half)
    ax.set_ylim(mid[1] - half, mid[1] + half)
    ax.set_zlim(mid[2] - half, mid[2] + half)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')


def plot_shape_only(ax, rv, tri, N_V):
    """Mode 1: plain shape."""
    polys = [[rv[t[0]], rv[t[1]], rv[t[2]]] for t in tri]
    pc = Poly3DCollection(polys, alpha=0.75,
                          facecolor='skyblue', edgecolor='k', linewidths=0.1)
    ax.add_collection3d(pc)
    ax.set_title('MeditRBC Shape')


def plot_shape_phi(ax, rv, tri, N_V, phi_dict):
    """Mode 2: shape coloured by phi."""
    polys = [[rv[t[0]], rv[t[1]], rv[t[2]]] for t in tri]
    fphi = face_scalar(tri, phi_dict, default=0.0)

    # Filter NaN
    valid = np.isfinite(fphi)
    if valid.any():
        vmin, vmax = fphi[valid].min(), fphi[valid].max()
    else:
        vmin, vmax = -1, 1
    if abs(vmax - vmin) < 1e-12:
        vmin -= 0.5; vmax += 0.5
    fphi = np.nan_to_num(fphi, nan=0.0)

    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.coolwarm
    fcolors = cmap(norm(fphi))

    pc = Poly3DCollection(polys, alpha=0.9)
    pc.set_facecolors(fcolors)
    pc.set_edgecolor('none')
    ax.add_collection3d(pc)

    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, shrink=0.6, label=r'$\varphi(r)$')
    ax.set_title(r'Concentration field $\varphi(r)$')


def plot_shape_phi_Q(ax, rv, tri, N_V, phi_dict, Q_dict, V_E, n_arrows=250):
    """Mode 3: shape coloured by phi, with director arrows from Q-tensor."""
    polys = [[rv[t[0]], rv[t[1]], rv[t[2]]] for t in tri]

    # Colour by phi if available, otherwise by S
    if phi_dict:
        fphi = face_scalar(tri, phi_dict, default=0.0)
        color_label = r'$\varphi(r)$'
    else:
        S_dict = {k: v[2] for k, v in Q_dict.items()}
        fphi = face_scalar(tri, S_dict, default=0.0)
        color_label = r'$S$ (order parameter)'

    valid = np.isfinite(fphi)
    if valid.any():
        vmin, vmax = fphi[valid].min(), fphi[valid].max()
    else:
        vmin, vmax = -1, 1
    if abs(vmax - vmin) < 1e-12:
        vmin -= 0.5; vmax += 0.5
    fphi = np.nan_to_num(fphi, nan=0.0)

    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.coolwarm
    fcolors = cmap(norm(fphi))

    pc = Poly3DCollection(polys, alpha=0.8)
    pc.set_facecolors(fcolors)
    pc.set_edgecolor('none')
    ax.add_collection3d(pc)

    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, shrink=0.6, label=color_label)

    # ----- Director arrows -----
    normals = compute_vertex_normals(rv, tri, N_V)
    T1, T2 = compute_tangent_frame(rv, tri, normals, V_E, N_V)

    # Sub-sample vertices for arrows
    vert_ids = sorted(Q_dict.keys())
    stride = max(1, len(vert_ids) // n_arrows)
    arrow_ids = vert_ids[::stride]

    # Determine a reasonable arrow length from mesh size
    pts = rv[1:N_V + 1]
    extent = (pts.max(axis=0) - pts.min(axis=0)).max()
    arrow_scale = extent / 30.0

    for vi in arrow_ids:
        q1, q2, S, D = Q_dict[vi]
        if S < 1e-6 or not np.isfinite(q1) or not np.isfinite(q2):
            continue
        alpha = 0.5 * np.arctan2(q2, q1)   # nematic director angle
        d = T1[vi] * np.cos(alpha) + T2[vi] * np.sin(alpha)
        dn = np.linalg.norm(d)
        if dn < 1e-15:
            continue
        d = d / dn * arrow_scale * min(S, 1.0)
        pos = rv[vi]
        ax.quiver(pos[0], pos[1], pos[2],
                  d[0], d[1], d[2],
                  color='black', linewidth=0.8, arrow_length_ratio=0.3)

    ax.set_title(r'Director $\mathbf{u}(r)$ + $\varphi(r)$')


# =============================================================
#  Main
# =============================================================

def main():
    parser = argparse.ArgumentParser(
        description='Visualize MeditRBC simulation output',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python visualize.py data/RMRC00.dat
  python visualize.py data/RMRC00.dat --phi data/Rfi00.dat
  python visualize.py data/RMRC00.dat --phi data/Rfi00.dat --Q data/RQ00.dat
  python visualize.py data/RMRC00.dat --Q data/RQ00.dat   (color by S)
""")
    parser.add_argument('shape', help='Shape file (Surface Evolver .dat)')
    parser.add_argument('--phi', help='Concentration field file')
    parser.add_argument('--Q', help='Q-tensor field file (q1 q2 S D)')
    parser.add_argument('-o', '--output', default='meditRBC_vis.png',
                        help='Output image file (default: meditRBC_vis.png)')
    parser.add_argument('--no-show', action='store_true',
                        help='Save image without opening display window')
    parser.add_argument('--arrows', type=int, default=250,
                        help='Approx. number of director arrows (default: 250)')
    parser.add_argument('--elev', type=float, default=25,
                        help='Elevation angle for 3D view (default: 25)')
    parser.add_argument('--azim', type=float, default=-60,
                        help='Azimuth angle for 3D view (default: -60)')
    args = parser.parse_args()

    if args.no_show:
        matplotlib.use('Agg')

    print(f'Reading shape: {args.shape}')
    rv_dict, V_E, E_F = read_shape(args.shape)
    rv, tri, N_V = build_mesh(rv_dict, V_E, E_F)
    print(f'  Vertices: {N_V}   Faces: {len(tri)}')

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev=args.elev, azim=args.azim)

    if args.phi is None and args.Q is None:
        # ---- Mode 1: shape only ----
        plot_shape_only(ax, rv, tri, N_V)

    elif args.Q is None:
        # ---- Mode 2: shape + phi ----
        print(f'Reading phi: {args.phi}')
        phi_dict = read_phi(args.phi)
        print(f'  phi vertices: {len(phi_dict)}')
        plot_shape_phi(ax, rv, tri, N_V, phi_dict)

    else:
        # ---- Mode 3: shape + (phi) + Q ----
        phi_dict = {}
        if args.phi:
            print(f'Reading phi: {args.phi}')
            phi_dict = read_phi(args.phi)
            print(f'  phi vertices: {len(phi_dict)}')
        print(f'Reading Q:   {args.Q}')
        Q_dict = read_Q(args.Q)
        print(f'  Q vertices: {len(Q_dict)}')
        plot_shape_phi_Q(ax, rv, tri, N_V, phi_dict, Q_dict, V_E,
                         n_arrows=args.arrows)

    set_equal_axes(ax, rv, N_V)
    plt.tight_layout()
    plt.savefig(args.output, dpi=150, bbox_inches='tight')
    print(f'Saved: {args.output}')

    if not args.no_show:
        plt.show()


if __name__ == '__main__':
    main()
