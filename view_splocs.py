import argparse
from itertools import count
import numpy as np
import h5py
from traits.api import HasTraits, Range, Instance, Bool, Int, on_trait_change
from traitsui.api import View, Item, HGroup, RangeEditor
from tvtk.api import tvtk
from tvtk.pyface.scene_editor import SceneEditor
from tvtk.common import configure_input, configure_input_data
from mayavi.tools.mlab_scene_model import MlabSceneModel
from mayavi.core.ui.mayavi_scene import MayaviScene
from pyface.timer.api import Timer

from util import veclen
from inout import load_splocs


class Visualization(HasTraits):
    component = Int(0)
    _max_component_index = Int()
    activation = Range(-1., 1.)
    oscillate = Bool(True)
    allow_negative = Bool(False)
    pd = Instance(tvtk.PolyData)
    normals = Instance(tvtk.PolyDataNormals)
    actor = Instance(tvtk.Actor)
    scene = Instance(MlabSceneModel, (), kw=dict(background=(1,1,1)))
    timer = Instance(Timer)

    def __init__(self, Xmean, tris, components):
        HasTraits.__init__(self)
        self._components = components
        self._max_component_index = len(components)
        self._Xmean = Xmean
        self.pd = tvtk.PolyData(points=Xmean, polys=tris)
        self.normals = tvtk.PolyDataNormals(splitting=False)
        configure_input_data(self.normals, self.pd)
        mapper = tvtk.PolyDataMapper(immediate_mode_rendering=True)
        self.actor = tvtk.Actor(mapper=mapper)
        configure_input(self.actor.mapper, self.normals)
        self.actor.mapper.lookup_table = tvtk.LookupTable(
            hue_range = (0.45, 0.6),
            saturation_range = (0., 0.8),
            value_range = (.6, 1.),
        )
        self.scene.add_actor(self.actor)
        self.timer = Timer(40, self.animate().next)

    def animate(self):
        for i in count():
            if self.oscillate:
                frame = i % 30
                alpha = np.sin(frame/30. * np.pi*2)
                if not self.allow_negative:
                    alpha = np.abs(alpha)
                self.activation = alpha
            yield

    @on_trait_change('activation, component')
    def update_plot(self):
        c = self._components[self.component]
        self.pd.points = self._Xmean + self.activation * c
        magnitude = veclen(c)
        self.pd.point_data.scalars = magnitude
        self.actor.mapper.scalar_range = (0, magnitude.max())
        self.scene.render()

    view = View(
        Item('scene', editor=SceneEditor(scene_class=MayaviScene),
             height=600, width=800, show_label=False),
        HGroup(
            Item('component', editor=RangeEditor(
                is_float=False, low=0, high_name='_max_component_index', mode='spinner')),
            'activation', 
            'oscillate', 
            'allow_negative',
        ),
        resizable=True, title="View SPLOC's",
    )

def main(component_hdf5_file):
    Xmean, tris, components, names = load_splocs(component_hdf5_file)

    visualization = Visualization(Xmean, tris, components)
    visualization.configure_traits()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Viewer for sparse localized deformation components')
    parser.add_argument('input_sploc_file')
    args = parser.parse_args()
    main(args.input_sploc_file)

