# Implementation of Modified Version of NTUPlace3 VLSI Cell Placement and Wirelength Minimization Algorithm
The repo should be cloned with `--recurse-submodules`. 
It is reccomended to run in a python venv, install the dependencies with `pip install -r requirements.txt`.

The entry point is `globalPlacement.py`. The input files are currently hard coded in the main function of this python file. Relavent LEF files for KSA16 test case are copied into the repo for development but eventually these should be pulled from OpenDB.
Running `globalPlacement.py` will run the placer and generate image files for all of the plots in the report.
