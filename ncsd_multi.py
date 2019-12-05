"""
==============
ncsd_multi.py
==============

A module to facilitate running ncsd with many nuclei.

It takes sets of inputs for multiple ncsd runs. Then, for each run:

1. It makes a folder so each run's output is contained, copies ``.exe`` file.
2. In each folder, creates a ``mfdp.dat`` file and a batch script.
3. Runs each bash script at the end, if desired. Otherwise you run them.

To run this file, manually change the values in this file then run::

    python ncsd_multi.py

"""
from os.path import realpath, dirname
from os import chdir
import sys
from multi_modules.data_structures import ManParams
from multi_modules.ncsd_multi_run import ncsd_multi_run
from multi_modules.data_checker import get_int_dir
from output_plotter import __file__ as output_plotter_path


ncsd_path = realpath("ncsd-it.exe")
working_dir = realpath("")
int_dir = realpath("../interactions/")
# you can also get int_dir from environment variable INT_DIR:
# int_dir = get_int_dir()
"""
Paths: change the values as needed, but don't remove the ``realpath``.
Empty string = current working directory, relative paths are relative to there.
"""

machine = "cedar"
"""Machine: must be one of "cedar", "summit", "local"."""

run = True
"""Do you want to run the code automatically? If so, set run=True"""

man_params = ManParams(
    # nucleus details:
    Z=3,  # number of protons
    N=5,  # number of neutrons
    hbar_omega=20,  # harmonic oscillator frequency
    N_1max=9,  # highest possible excited state of 1 nucleon
    N_12max=10,  # highest possible state of 2 nucleons, added
    # if 3-body, change this and interaction name, otherwise leave them be
    N_123max=11,  # highest possible state of 3 nucleons, added

    # interaction names, these files must be within int_dir
    two_body_interaction="some_tbme_file",
    three_body_interaction="some_three_body_file",
    # potential is just for naming purposes, does not affect calculations
    potential_name="some_potential",

    # computation-related parameters:
    Nmax_min=0,  # Nmax_min and Nmax_max control how long the program runs
    Nmax_max=10,  # e.g. 0 - 8 gives eigenvectors for Nmax = 0,2,...,8
    Nmax_IT=6,  # Nmax for Importance Truncation
    interaction_type=2,  # make sure abs(interaction_type)==3 for 3-body
    n_states=1,  # number of final states (= number of energy values)
    iterations_required=10,  # number of iterations required in Lanczos step
    irest=0,  # restart? 4 = yes, 0 = no
    nhw_restart=-1,  # not sure what this one does
    kappa_points=4,  # number of kappa values
    # make sure you keep this next one written as a string!
    kappa_vals="2.0 3.0 5.0 10.0",  # values for kappa, in increasing order
    kappa_restart=-1,  # -1 for false, else some value between 1 and 4
    saved_pivot="F",  # "T" or "F", whether or not to use the saved pivot

    # machine-related parameters:
    time="0 8 0",  # runtime, "days hours minutes". Will be formatted later
    mem=80.0,  # memory, in GB
    n_nodes=1024  # number of nodes
)
"""
Manual Parameters:

- specify all as single parameter or list []
- if you use a list, we'll make N ncsd runs, where N is the maximum length
  of a list made here.
- list of attributes is in multi_modules/data_structures.py, ``man_keys``.
- default parameters are at the bottom of multi_modules/data_structures.py
"""


def run(man_params, int_dir, ncsd_path, working_dir, machine, run=True):
    """
    Runs NCSD code. The real heavy lifting is done in ncsd_multi_run.py.

    man_params:
        a data_structures.ManParams instance,
        containing all the manual parameters you want to enter

    int_dir:
        a path to the directory where you keep interaction files

    ncsd_path:
        the path to the ncsd executable file

    working_dir:
        the path to the directory in which you want to create
        new directories for each run

    machine:
        the name of the machine you're using, e.g. "cedar", "summit"

    Returns:
        Nothing, but does create files and optinally runs batch scripts.
    """
    sys.path.append(output_plotter_path)
    # sys.tracebacklimit = 0  # Suppresses tracebacks, if you want that
    paths = [int_dir, ncsd_path, working_dir, output_plotter_path]
    # set run=True to run all batch scripts
    ncsd_multi_run(man_params, paths, machine, run=run)


if __name__ == "__main__":
    run(man_params, int_dir, ncsd_path, working_dir, machine, run=run)
