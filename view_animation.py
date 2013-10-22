import argparse
import h5py
from itertools import count
from mayavi import mlab
from tvtk.api import tvtk

def main(hdf5_animation_file):
    weights = None
    with h5py.File(hdf5_animation_file, 'r') as f:
        verts = f['verts'].value
        tris = f['tris'].value
        if 'weights' in f:
            weights = f['weights'].value

    pd = tvtk.PolyData(points=verts[0], polys=tris)
    normals = tvtk.PolyDataNormals(input=pd, splitting=False)
    actor = tvtk.Actor(mapper=tvtk.PolyDataMapper(input=normals.output))
    actor.property.set(edge_color=(0.5, 0.5, 0.5), ambient=0.0,
                       specular=0.15, specular_power=128., shading=True, diffuse=0.8)

    fig = mlab.figure(bgcolor=(1,1,1))
    fig.scene.add_actor(actor)

    @mlab.animate(delay=40, ui=False)
    def animation():
        for i in count():
            if weights is not None:
                w_str = ",".join(["%0.2f"] * weights.shape[1])
                print ("Frame %d Weights = " + w_str) % tuple([i] + weights[i].tolist())
            frame = i % len(verts)
            pd.points = verts[frame]
            fig.scene.render()
            yield

    a = animation()
    fig.scene.z_minus_view()
    mlab.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Animation viewer for hdf5 mesh animationfiles '
                    '(use import scripts to convert other formats to hdf5)')
    parser.add_argument('input_filename')
    args = parser.parse_args()
    main(args.input_filename)

