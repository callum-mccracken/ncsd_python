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
from os.path import realpath
import sys
from multi_modules.data_structures import ManParams
from multi_modules.ncsd_multi_run import ncsd_multi_run
from multi_modules.data_checker import get_int_dir, get_work_dir
from output_plotter import __file__ as output_plotter_path


ncsd_path = realpath("ncsd-it.exe")
working_dir = realpath("")
int_dir = realpath("../interactions/")
# you can also get these directories from
# the environment variables INT_DIR and WORK_DIR:
# int_dir = get_int_dir()
# working_dir = get_work_dir()
"""
Paths: change the values as needed, but don't remove the ``realpath``
(it converts relative --> absolute paths).
Relative paths = relative to current working directory.
"""

machine = "cedar"
"""machine must be one of "cedar", "summit", "local"."""

run = True
"""Do you want to run jobs automatically? If so, set run=True"""

man_params = ManParams(
    Z=3,  # number of protons
    N=5,  # number of neutrons
    hbar_omega=20,  # harmonic oscillator frequency
    N_1max=9,  # highest number of harmonic oscillator quanta for 1 nucleon
    N_12max=10,  # highest number of harmonic oscillator quanta for 2 nucleons
    # Nmax, in contrast to N_1max / N_12max, is the max number of oscillator
    # quanta allowed above the lowest Pauli-allowed state for the nucleus
    Nmax_min=0,  # Nmax_min and Nmax_max are lower/upper bounds for Nmax
    Nmax_max=10,  # e.g. 0 - 8 gives eigenvectors for Nmax = 0,2,...,8
    Nmax_IT=6,  # Nmax value at which importance truncation starts
    n_states=1,  # number of final states (= number of energy values)
    iterations_required=10,  # number of iterations required in Lanczos step
    irest=0,  # restart? 4 = yes, 0 = no
    nhw_restart=-1,  # nhw restart? -1 or 1, I think
    kappa_points=4,  # number of kappa values
    # note kappa values are x10^-4
    kappa_vals="2.0 3.0 5.0 10.0",  # kappa values, increasing, use a string!
    kappa_restart=-1,  # -1 for false, else some value between 1 and 4
    saved_pivot="F",  # "T" or "F", whether or not to use the saved pivot
    time="0 8 0",  # runtime, "days hours minutes". Will be formatted later
    mem=80.0,  # memory, in GB
    n_nodes=1024,  # number of nodes
    potential_name="n3lo-NN3Nlnl-srg2.0",  # identifier for output files

    # 2-body interaction filename, must be within int_dir (defined above)
    two_body_interaction="some_tbme_file",

    interaction_type=2,  # 2: 2-body, 3: 3-body, -3: 3-body with 2-body IT
    # if 3-body interaction, change these two, if not they will be ignored
    N_123max=11,  # highest number of harmonic oscillator quanta for 3 nucleons
    three_body_interaction="some_three_body_file",  # must be within int_dir!
)
"""
Manual Parameters:

- specify all as single parameter or list []
- if you use a list, we'll make N ncsd runs, where N is the maximum length
  of a list made here.
- list of attributes is in multi_modules/data_structures.py, ``man_keys``.
- default parameters are at the bottom of multi_modules/data_structures.py
"""


def run_ncsd(man_params, int_dir, ncsd_path, working_dir, machine, run=True):
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
    run_ncsd(man_params, int_dir, ncsd_path, working_dir, machine, run=run)
