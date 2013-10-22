import sys
import h5py
import pylab as pl

if __name__ == '__main__':
    with h5py.File(sys.argv[1], 'r') as f:
        weights = f['weights'].value

    pl.figure(figsize=(4, 12))
    for i in xrange(weights.shape[1]):
        pl.subplot(weights.shape[1], 1, i+1)
        pl.plot(weights[:,i])
        pl.xticks([])
        pl.yticks([])
    pl.show()

