# ncsd_python

Python modules to make your life easy while running `ncsd-it.exe`!

Main script is `ncsd_multi.py`. You feed it a list of parameters,
and it'll run multiple ncsd jobs for you, then plot the output.

## Getting Started

`git clone` this, or if you're on Cedar you can find a clone in `exch`.

To get the latest version, do:

`git pull origin master`

or for a hard reset, do

`git fetch --all`

`git reset --hard origin/master`

Then you'll need to open `ncsd_multi.py` and enter a bunch of parameters:

- `ncsd_path = realpath("ncsd-it.exe")`
  - edit this so it's the path to your ncsd executable file
  - relative paths are okay inside the `realpath()` function
- `working_dir = realpath("")`
  - make this the path to the directory where to want to complete the ncsd runs
- `int_dir = realpath("../interactions")`
  - make this the path to your interactions directory

```
    manual_params = ManParams(
        ...
        Z = 3,  # number of protons
        N = [5,6],  # number of neutrons
        hbar_omega = 20,  # harmonic oscillator frequency
        N_1max = 9,  # highest possible excited state of 1 nucleon
        N_12max = 10,  # highest possible state of 2 nucleons, added
        ...
    )
```
  - enter your calculation-related parameters here

- `ncsd_multi_run(man_params, paths, machine run=True)`
  - set `run` to `False` if you don't want to run all batch files

- There are other parameters which are set by default, e.g. `iclmb`
  - Those can be changed by going to the very bottom of `data_structures.py`
    and editing the inputs to `DefaultParamsObj`

### To do multiple runs, simply specify at least one variable as a list.
e.g. `Z = [3,4,5]` will make 3 different runs.

Make sure that your variables are well ordered if more than one is changing,
runs are created by list index!

If you have variables that change, but less than others,
keep in mind that lists are extended by copying their last entry.

So if you have

`Z = [1,2]`
`N = [1,2,3]`
`Nhw = 1`

you'll get 3 runs, with parameters

`Z = [1,2,2]`
`N = [1,2,3]`
`Nhw = [1,1,1]`

Note: make sure to edit the 3-body parameters if `abs(interaction_type) == 3`.

## Output Plotting

Once each job is run, `output_plotter.py` will be called, which plots the output
from ncsd in `xmgrace` format by default (but can produce other output).

`output_plotter.py` can be run on its own too, open the file and have a look
at the parameters there for extra details.

## Prerequisites

- Python (3.7.4 ideally, other versions may work)
- `numpy`, install with `pip install --user numpy`
- `matplotlib`, install with `pip install --user matplotlib`
  
  
- NCSD executable
- interaction files
- `grace_spectra_plotter.exe`
- If you want to get data from BNL:
  - `tabula`, install with `pip install --user tabula-py`
  - this is recommended over TUNL, but both are pretty experimental
- If you want to get data from TUNL:
  - Java SDK
    - if you type `java` in the terminal and get sensible output, you probably have it
  - A few python libraries:
    - `tabula`, install with `pip install --user tabula-py`
    - `requests`, install with `pip install --user requests`
    - `lxml`, install with `pip install --user lxml`
    - `matplotlib`, install with `pip install --user matplotlib`
