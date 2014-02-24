import sys
import h5py
import pylab as pl

if __name__ == '__main__':
    weights = []
    for filename in sys.argv[1:]:
        with h5py.File(filename, 'r') as f:
            weights.append(f['weights'].value)

    pl.figure(figsize=(4, 12))
    for i in xrange(weights[0].shape[1]):
        pl.subplot(weights[0].shape[1], 1, i+1)
        for w in weights:
            pl.plot(w[:,i])
        pl.xticks([])
        pl.yticks([])
    pl.show()

