=======
picasso
=======

Python ICe AnalySis SOftware!

Python package for analyzing GRACE data.
TODO: detailed description.


Environment
===========

To set up the environment install Miniconda and use:

::

  conda env create -f environment.yml
  conda activate picasso
  python setup.py develop

Note: change "develop" to "install" if you are not planning on changing the source code.

Open PICASSO GUI
================

The software can be accessed via a GUI. In a terminal window run

::

  picasso-gui

A less complex command line interface might follow in the future.

Modify the GUI
==============

The GUI was created with QT Designer. Appropriate changes can be made in the picasso_design.ui file located at /picasso/src/data/gui_data. 
The modified ui-file needs to be converted to a Python module using the script ui2py.sh located in the same directory.

Note
====

This project has been set up using PyScaffold 3.1. For details and usage
information on PyScaffold see https://pyscaffold.org/.


