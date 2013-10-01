splocs
======

A reference implementation of "Sparse Localized Deformation Components". Notice that this reference implementation does not include all the features from the paper, e.g. it currently misses user priors and convergence checks. The implementation is meant for documentation purposes. The aim is to keep the code clean, short, and well documented. The code is still able to run, of course ;-) In this way, I hope that future researchers will be able to quickly come up with new applications of splocs and I encourage everyone to try variations of splocs for their own research or even commercial problems. 

The code currently doesn't come with a dataset, but there are scripts to import OFF and PLY mesh sequences easily. 

The main algorithm presented in the paper is in "splocs.py"

I would be very happy to hear about cool applications built with this code. Any feedback is also very appreciated.

Requirements
------------

 - Python
 - numpy
 - scipy

For visualizing the components, additionally
 - mayavi2

Usage Example # 1 with already preprocessed data
 - Look at animation
        >>> python view_animation.py data/face_lizhang_5000.h5
 - Find sparse components (takes some time)
        >>> python sploc.py data/face_lizhang_5000.h5 /tmp/lizhang_splocs.h5
        Notice that right now, you can only change the algorithm parameters in the file sploc.py (around line 30)
 - View sparse components
        >>> python view_splocs.py /tmp/lizhang_splocs.h5

The sploc algorithm requires the input files pre-processed and aligned, for own data, use the following process:
 - Convert and preprocess the mesh sequence (OFF and PLY import script available):
        >>> python import_off.py  "~/data/volker/*.off" /tmp/volker.h5
        >>> python align_rigid.py /tmp/volker.h5 /tmp/volker_aligned.h5 -r first
 - Look at animation
        >>> python view_animation.py /tmp/volker_aligned.h5

Licensing
---------

This code is licenced under the MIT License. If you use this code in publications, please cite:

        Thomas Neumann, Kiran Varanasi, Stephan Wenger, Markus Wacker, Marcus Magnor, and Christian Theobalt
        Sparse Localized Deformation Components
        ACM Transactions on Graphics 32 (6), 2013 (Proceedings of SIGGRAPH Asia)

Here is the Bibtex entry:
@article{neumann2013splocs,
    author = {Neumann, Thomas and Varanasi, Kiran and Wenger, Stephan and Wacker, Markus and Magnor, Marcus and Theobalt, Christian},
    title = {Sparse Localized Deformation Components},
    volume = {32},
    number = {6},
    journal = {{ACM} Transactions on Graphics (Proc. of Siggraph Asia)},
    month = nov,
    year = {2013}
}

