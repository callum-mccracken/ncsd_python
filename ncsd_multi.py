"""
ncsd_multi.py: facilitates running ncsd-it.exe with many nuclei

- takes sets of inputs for ncsd-it.exe runs
- then, for each run:
    - makes a folder so each run's output is contained, copies exe file
    - makes a bunch of mfdp.dat files to complete each run
    - runs each process at the end, if desired
- more details in README.md
"""
from os.path import realpath, dirname
from os import chdir
import sys
from multi_modules.data_structures import ManParams
from multi_modules.ncsd_multi_run import ncsd_multi_run
from multi_modules.data_checker import get_int_dir
from output_exporter import __file__ as output_plotter_path

# sys.tracebacklimit = 0  # If debugging comment this out! Suppresses tracebacks

# change these to suit your needs, but don't remove the "realpath"
# empty string = current working directory, relative paths = relative to cwd
ncsd_path = realpath("ncsd-it.exe")
working_dir = realpath("")
int_dir = realpath("../interactions/")
# you can also get int_dir from environment variable INT_DIR:
# int_dir = get_int_dir()

# set machine name, make sure it's valid
machine = "local"
assert machine in ["cedar", "summit", "local"]

# PARAMETERS -- specify all as single parameter or list []
man_params = ManParams(
    # nucleus details:
    Z=1,  # number of protons
    N=1,  # number of neutrons
    hbar_omega=20,  # harmonic oscillator frequency
    N_1max=1,  # highest possible excited state of 1 nucleon
    N_12max=2,  # highest possible state of 2 nucleons, added
    # if 3-body, change this and interaction name, otherwise just leave them be
    N_123max=11,  # highest possible state of 3 nucleons, added

    # interaction names, these files must be within int_dir
    two_body_interaction="TBME.int",
    three_body_interaction="none",
    # potential is just for naming purposes, does not affect calculations
    potential_name="test_potential",

    # computation-related parameters:
    Nmax_min=0,  # Nmax_min and Nmax_max control how long the program runs
    Nmax_max=10,  # e.g. 0 - 8 will give you eigenvectors for Nmax = 0,2,...,8
    Nmax_IT=6,  # Nmax for Importance Truncation
    interaction_type=2,  # make sure abs(interaction_type)==3 for 3-body
    n_states=1,  # number of final states (= number of energy values)
    iterations_required=10,  # number of iterations required in Lanczos step
    irest=0,  # restart? 1 = yes, 0 = no
    nhw_restart=-1,  # not sure what this one does
    kappa_points=4,  # number of kappa values
    # make sure you keep this next one written as a string!
    kappa_vals="2.0 3.0 5.0 10.0",  # values for kappa, in increasing order
    kappa_restart=-1,  # -1 for false, else some value between 1 and 4
    saved_pivot="F",  # "T" or "F", whether or not to use the saved pivot

    # machine-related parameters:
    time="0 8 0",  # runtime, "days hours minutes". It will be formatted later
    mem=1.0,  # memory, in GB
    n_nodes=1024  # number of nodes
)

# default parameters can be found at the bottom of data_structures.py
# (which is in the sub_modules directory)

paths = [int_dir, ncsd_path, working_dir, output_plotter_path]

# set run=True to run all batch scripts
ncsd_multi_run(man_params, paths, machine, run=True)
