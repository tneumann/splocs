splocs
======

A reference implementation of "Sparse Localized Deformation Components", a paper accepted for SIGGRAPH Asia 2013.  More information on the paper, including the videos, can be found on the [project page](http://www.drematrix.de/?portfolio=english-sparse-localized-deformation-components&lang=en).

The main algorithm presented in the paper is in "sploc.py"

I would be very happy to hear about cool applications built with this code, and any feedback is very appreciated.

Requirements
------------

To run the code, you will need
 - Python
 - numpy
 - scipy

For visualizing the components and animations, additionally
 - mayavi2

In Ubuntu, these can be easily installed by

        >>> apt-get install python python-numpy python-scipy mayavi2

Usage
-----

The sploc algorithm requires the input files pre-processed and aligned, for own data, use the following process:

 - Convert and preprocess the mesh sequence (OFF and PLY import script available):

        >>> python import_off.py  "~/data/volker/*.off" /tmp/volker.h5

        >>> python align_rigid.py /tmp/volker.h5 /tmp/volker_aligned.h5 -r first

 - Look at animation

        >>> python view_animation.py /tmp/volker_aligned.h5

Here is an example of processing the data

 - Find sparse components (takes some time)

        >>> python sploc.py /tmp/volker_aligned.h5 /tmp/volker_splocs.h5

    Notice that right now, you can only change the algorithm parameters in the file sploc.py (around line 30)

 - View sparse components

        >>> python view_splocs.py /tmp/volker_splocs.h5


Obtaining the datasets used in the paper
----------------------------------------

Here is how you can obtain the datasets that are used in the paper:

**[Zhang et al.2004]** Facial performance capture data. You can get this here, as a PLY sequence: http://grail.cs.washington.edu/software-data/stfaces/index.html

**[Valgaerts et al. 2012]** High-resolution facial performance capture data. You can request the data here: http://gvv.mpi-inf.mpg.de/projects/FaceCap/

**[Beeler et al. 2011]** High-resolution facial performance capture data. The instructions how to get it are here: http://graphics.ethz.ch/publications/papers/paperBee11.php

**[Hasler et al. 2009]** Full body scans in correspondence. Request here: http://gvvperfcapeva.mpi-inf.mpg.de/public/ScanDB/

**[Neumann et al. 2013]** Captured shoulder and arm muscles. Please ask by mail, see the project page http://www.drematrix.de/publications/capture-and-statistical-modeling-of-arm-muscle-def/ This data is not in OFF/PLY Format and needs a separate pose normalization step that is not part of this source code distribution.

Limitations
-----------

Notice that this reference implementation does not include all the features from the paper, e.g. it currently misses user priors and convergence checks. The implementation is meant for documentation purposes. The aim is to keep the code clean, short, and well documented. The code is still able to run, of course ;-) In this way, I hope that future researchers will be able to quickly come up with new applications of splocs and I encourage everyone to try variations of splocs for their own research or even commercial problems. 

The code currently doesn't come with a dataset, but there are scripts to import OFF and PLY mesh sequences easily. 

License
-------

This code is licenced under the MIT License. If you use this code in publications, please cite:

        Thomas Neumann, Kiran Varanasi, Stephan Wenger, Markus Wacker, Marcus Magnor, and Christian Theobalt
        Sparse Localized Deformation Components
        ACM Transactions on Graphics 32 (6), 2013 (Proceedings of SIGGRAPH Asia)

Bibtex entry:
```
@article{neumann2013splocs,
    author = {Neumann, Thomas and Varanasi, Kiran and Wenger, Stephan and Wacker, Markus and Magnor, Marcus and Theobalt, Christian},
    title = {Sparse Localized Deformation Components},
    volume = {32},
    number = {6},
    journal = {{ACM} Transactions on Graphics (Proc. of Siggraph Asia)},
    month = nov,
    year = {2013}
}
```


[![Bitdeli Badge](https://d2weczhvl823v0.cloudfront.net/tneumann/splocs/trend.png)](https://bitdeli.com/free "Bitdeli Badge")

