from os import path
from glob import glob
from cStringIO import StringIO
import numpy as np
import h5py
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix

from util import sort_nicely, veclen, filter_reindex


def convert_sequence_to_hdf5(filename_pattern, loader_function, hdf_output_file):
    verts_all = []
    tris = None
    files = glob(path.expanduser(filename_pattern))
    sort_nicely(files)
    for i, f in enumerate(files):
        print "loading file %d/%d [%s]" % (i+1, len(files), f)
        verts, new_tris = loader_function(f)
        if tris is not None and new_tris.shape != tris.shape and new_tris != tris:
            raise ValueError, "inconsistent topology between meshes of different frames"
        tris = new_tris
        verts_all.append(verts)

    verts_all = np.array(verts_all, np.float32)
    verts_all, tris, _, verts_mean, verts_scale = preprocess_mesh_animation(verts_all, tris)

    with h5py.File(hdf_output_file, 'w') as f:
        f.create_dataset('verts', data=verts_all, compression='gzip')
        f['tris'] = tris
        f.attrs['mean'] = verts_mean
        f.attrs['scale'] = verts_scale

    print "saved as %s" % hdf_output_file

def preprocess_mesh_animation(verts, tris):
    """ 
    Preprocess the mesh animation:
        - removes zero-area triangles
        - keep only the biggest connected component in the mesh
        - normalize animation into -0.5 ... 0.5 cube
    """
    print "Vertices: ", verts.shape
    print "Triangles: ", verts.shape
    assert verts.ndim == 3
    assert tris.ndim == 2
    # check for zero-area triangles and filter
    e1 = verts[0, tris[:,1]] - verts[0, tris[:,0]]
    e2 = verts[0, tris[:,2]] - verts[0, tris[:,0]]
    n = np.cross(e1, e2)
    tris = tris[veclen(n) > 1.e-8]
    # remove unconnected vertices
    ij = np.r_[np.c_[tris[:,0], tris[:,1]], 
               np.c_[tris[:,0], tris[:,2]], 
               np.c_[tris[:,1], tris[:,2]]]
    G = csr_matrix((np.ones(len(ij)), ij.T), shape=(verts.shape[1], verts.shape[1]))
    n_components, labels = connected_components(G, directed=False)
    if n_components > 1:
        size_components = np.bincount(labels)
        if len(size_components) > 1:
            print "[warning] found %d connected components in the mesh, keeping only the biggest one" % n_components
            print "component sizes: "
            print size_components
        keep_vert = labels == size_components.argmax()
    else:
        keep_vert = np.ones(verts.shape[1], np.bool)
    verts = verts[:, keep_vert, :]
    tris = filter_reindex(keep_vert, tris[keep_vert[tris].all(axis=1)])
    # normalize triangles to -0.5...0.5 cube
    verts_mean = verts.mean(axis=0).mean(axis=0)
    verts -= verts_mean
    verts_scale = np.abs(verts.ptp(axis=1)).max()
    verts /= verts_scale
    print "after preprocessing:"
    print "Vertices: ", verts.shape
    print "Triangles: ", verts.shape
    return verts, tris, ~keep_vert, verts_mean, verts_scale

def load_ply(filename):
    try:
        from enthought.tvtk.api import tvtk
    except ImportError:
        try:
            from tvtk.api import tvtk
        except ImportError:
            print "Reading PLY files requires TVTK. The easiest way is to install mayavi2"
            print "(e.g. on Ubuntu: apt-get install mayavi2)"
            raise
    reader = tvtk.PLYReader(file_name=filename)
    reader.update()
    polys = reader.output.polys.to_array().reshape((-1, 4))
    assert np.all(polys[:,0] == 3)
    return reader.output.points.to_array(), polys[:,1:]

def load_off(filename, no_colors=False):
    lines = open(filename).readlines()
    lines = [line for line in lines if line.strip() != '' and line[0] != '#']
    assert lines[0].strip() in ['OFF', 'COFF'], 'OFF header missing'
    has_colors = lines[0].strip() == 'COFF'
    n_verts, n_faces, _ = map(int, lines[1].split())
    vertex_data = np.loadtxt(
        StringIO(''.join(lines[2:2 + n_verts])), 
        dtype=np.float)
    if n_faces > 0:
        faces = np.loadtxt(StringIO(''.join(lines[2+n_verts:])), dtype=np.int)[:,1:]
    else:
        faces = None
    if has_colors:
        colors = vertex_data[:,3:].astype(np.uint8)
        vertex_data = vertex_data[:,:3]
    else:
        colors = None
    if no_colors:
        return vertex_data, faces
    else:
        return vertex_data, colors, faces

def save_off(filename, vertices=None, faces=None):
    if vertices is None:
        vertices = []
    if faces is None:
        faces = []
    with open(filename, 'w') as f:
        f.write("OFF\n%d %d 0\n" % (len(vertices), len(faces)))
        if len(vertices) > 1:
            np.savetxt(f, vertices, fmt="%f %f %f")
        if len(faces) > 1:
            for face in faces:
                fmt = " ".join(["%d"] * (len(face) + 1)) + "\n"
                f.write(fmt % ((len(face),) + tuple(map(int, face))))

def load_splocs(component_hdf5_file):
    with h5py.File(component_hdf5_file, 'r') as f:
        tris = f['tris'].value
        Xmean = f['default'].value
        names = sorted(list(set(f.keys()) - set(['tris', 'default'])))
        components = np.array([
            f[name].value - Xmean 
            for name in names])
    return Xmean, tris, components, names
