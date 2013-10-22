import numpy as np
import h5py
import argparse

from inout import load_splocs


def main(input_animation_file, input_sploc_file, output_animation_file):
    with h5py.File(input_animation_file, 'r') as f:
        verts = f['verts'].value
        tris = f['tris'].value
    Xmean, tris_splocs, C, names = load_splocs(input_sploc_file)

    assert tris.shape == tris_splocs.shape, "animation should have same topology"
    assert np.all(tris == tris_splocs), "animation should have same topology"

    # fit weights
    X = (verts - Xmean)
    # pre-invert component subspace so that all frames can be optimized very efficiently
    Cflat = C.reshape(C.shape[0], -1)
    Cinv = np.linalg.pinv(np.dot(Cflat, Cflat.T) + np.eye(Cflat.shape[0]) * 1.e-5)
    # optimize for weights of each frame
    W = np.array(
        [np.dot(Cinv, np.dot(Cflat, x.ravel())) for x in X])
    # reconstruct animation
    verts_reconstructed = np.tensordot(W, C, (1, 0)) + Xmean[np.newaxis]
    # save
    with h5py.File(output_animation_file, 'w') as f:
        f['verts'] = verts_reconstructed
        f['tris'] = tris
        f['weights'] = W


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Find weights for given input animation and given sparse components, output the reconstructed animation including the weights. Performs a simple least-squares fit of the weights and thus may yield unconstrained weights (e.g. not in the usual range between 0..1).')
    parser.add_argument('input_animation_file')
    parser.add_argument('input_sploc_file')
    parser.add_argument('output_animation_file',
                        help='Output fitted animation file (will also save the component weights)')
    args = parser.parse_args()
    main(args.input_animation_file, args.input_sploc_file, args.output_animation_file)
