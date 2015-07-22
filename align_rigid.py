import argparse
import numpy as np
import h5py

from util import transform


def find_rbm_procrustes(frompts, topts):
    """
    Finds a rigid body transformation M that moves points in frompts to the points in topts
    that is, it finds a rigid body motion [ R | t ] with R \in SO(3)

    This algorithm first approximates the rotation by solving
    the orthogonal procrustes problem.
    """
    # center data
    t0 = frompts.mean(0)
    t1 = topts.mean(0)
    frompts_local = frompts - t0
    topts_local = topts - t1
    # find best rotation - procrustes problem
    M = np.dot(topts_local.T, frompts_local)
    U, s, Vt = np.linalg.svd(M)
    R = np.dot(U, Vt)
    if np.linalg.det(R) < 0:
        R *= -1
    T0 = np.eye(4)
    T0[:3,:3] = R
    T0[:3, 3] = t1 - np.dot(R, t0)
    return T0


def test_find_rbm_procrustes():
    for it in xrange(100):
        R = np.linalg.qr(np.random.uniform(-1, 1, size=(3, 3)))[0]
        if np.linalg.det(R) < 0:
            R *= -1
        t = np.random.uniform(-2, 2, size=3)
        M = np.eye(4)
        M[:3, :3] = R
        M[:3,  3] = t
        N = np.random.randint(3, 1000)
        frompts = np.random.random((N, 3))
        topts = transform(frompts, M)
        M_pro = find_rbm_procrustes(frompts, topts)
        np.testing.assert_almost_equal(M_pro, M)


def main(input_hdf5_file, output_hdf5_file):
    data = h5py.File(input_hdf5_file, 'r')
    verts = data['verts'].value
    tris = data['tris'].value

    v0 = verts[0]
    verts_new = []
    for i, v in enumerate(verts):
        print "frame %d/%d" % (i+1, len(verts))
        M = find_rbm_procrustes(v, v0)
        verts_new.append(transform(v, M))
    verts = np.array(verts_new, np.float32)

    with h5py.File(output_hdf5_file, 'w') as f:
        f.create_dataset('verts', data=verts, compression='gzip')
        f['tris'] = tris

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Align mesh animation rigidly')
    parser.add_argument('input_animation_file')
    parser.add_argument('output_animation_file')
    args = parser.parse_args()
    main(args.input_animation_file, args.output_animation_file)

