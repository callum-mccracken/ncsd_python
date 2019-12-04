"""
Once we grab input from ncsd_multi.py, this module actually processes it.
"""
from os import system, chdir, mkdir, symlink
from os.path import realpath, join, exists, relpath, dirname
from shutil import rmtree

from .data_structures import ManParams
from .parameter_calculations import calc_params, nucleus
from .data_checker import manual_input_check
from .file_manager import MFDP, CedarBatch, SummitBatch, LocalBatch, Defaults


def prepare_input(m_params):  # m_params for manual params
    """
    Takes manual input m_params (a ManParams instance),
    and parses it into a format that is usable by the functions
    that make the run directories.

    Makes lists of length num_runs (longest length list in m_params) with
    duplicates of the last element if there aren't enough.

    For example, if num_runs is 4...

    1 --> [1,1,1,1]

    [1,2] --> [1,2,2,2]

    This may be useful to know, in case if you thought

    [1,2] --> [1,1,2,2] (wrong!)

    """
    print("preparing input to be written to files")
    m_dict = m_params.param_dict()
    # make one set of inputs for each run

    # make all parameters into lists
    for key, value in m_dict.items():
        if type(value) != list:
            m_dict[key] = [value]

    #  now length of longest list = number of unique runs
    num_runs = len(max(list(m_dict.values()), key=len))

    for key, value in m_dict.items():
        if len(value) != num_runs:  # it's also guaranteed < num_runs
            # append the last element until it's the right length
            for i in range(len(value), num_runs):
                m_dict[key].append(value[-1])

    # now make one dict for each mfdp file:
    dict_list = []
    for i in range(num_runs):
        new_dict = {}
        for key, value in m_dict.items():
            new_dict[key] = value[i]
        dict_list.append(new_dict)
    return dict_list


def create_dirs(defaults, dict_list, paths, machine):
    """
    Creates directories which hold ncsd run files.

    Returns paths to batch scripts in each directory. Run these to run ncsd.

    defaults:
        a data_structures.DefaultParams object, which stores default values

    dict_list:
        a list of ManParams objects, created from the original ManParams object

    paths:
        a list containing
        [interactions_dir, ncsd_path, working_dir, output_plotter_path]

    machine:
        name of machine being used for this calculation, e.g. "cedar"
    """
    print("creating directories to store run files")

    def populate_dir(defaults, man_params, paths, machine):
        """
        Creates directories with mfdp.dat and batch files,
        using both manual input and default parameters.

        defaults:
            a data_structures.DefaultParams object, which stores default values

        man_params:
            a ManParams object containing variables for this ncsd run

        paths:
            a list containing
            [interactions_dir, ncsd_path, working_dir, output_plotter_path]

        machine:
            name of machine being used for this calculation, e.g. "cedar"
        """
        ncsd_path = paths[1]
        working_dir = paths[2]

        chdir(working_dir)

        # make a directory for run
        run_name = nucleus(man_params.Z, man_params.N)
        run_dir = realpath(join(working_dir, run_name))
        # ensure we don't overwrite
        if exists(run_dir):
            new_name = input(
                "Run '"+run_name+"' already exists. \n"
                "Enter new name, or hit enter to overwrite: ")
            if new_name:
                run_dir = realpath(join(working_dir, new_name))
            else:
                #  remove it and start from scratch
                rmtree(run_dir)
        print("making run directory "+run_dir)
        mkdir(run_dir)

        # now actually calculate the parameters to write out
        [mfdp_params, batch_params] = calc_params(
            run_dir, paths, man_params, defaults.params, machine)

        # be sure that all the batch files actually know where their exe is
        batch_params.ncsd_path = realpath(join(run_dir, "ncsd-it.exe"))

        print("writing files")
        # copy ncsd-it.exe
        symlink(ncsd_path, batch_params.ncsd_path)

        # convert interaction files to relative paths too
        mfdp_params.two_body_interaction = relpath(
            mfdp_params.two_body_interaction, run_dir)
        mfdp_params.three_body_interaction = relpath(
            mfdp_params.three_body_interaction, run_dir)

        # write mfdp.dat file
        mfdp_path = realpath(join(run_dir, "mfdp.dat"))
        MFDP(filename=mfdp_path, params=mfdp_params).write()

        # before writing bacth file, convert ncsd_path to relative path
        batch_params.ncsd_path = relpath(batch_params.ncsd_path, run_dir)
        # ensure ncsd path has a ./ if needed
        if "/" not in batch_params.ncsd_path:
            batch_params.ncsd_path = "./" + batch_params.ncsd_path

        # write batch file
        batch_path = realpath(join(run_dir, "batch_ncsd"))
        # but first make sure to convert output_plotter to a relative path
        batch_params.output_plotter = relpath(
            batch_params.output_plotter, run_dir)

        if machine == "local":
            LocalBatch(filename=batch_path, params=batch_params).write()
        elif machine == "cedar":
            CedarBatch(filename=batch_path, params=batch_params).write()
        elif machine == "summit":
            SummitBatch(filename=batch_path, params=batch_params).write()
        else:
            raise ValueError("Invalid machine!")

        # then tell the program where it is so we can run it later
        return batch_path

    # for each set of inputs
    batch_paths = []
    for man_params_dict in dict_list:
        man_params = ManParams(**man_params_dict)
        batch_paths.append(
            populate_dir(defaults, man_params, paths, machine)
        )
    # return list of paths to be run
    return batch_paths


def ncsd_multi_run(man_params, paths, machine, run=True):
    """
    Run ncsd multiple times with given parameters.

    man_params:
        a ManParams object with possible lists of values for all ncsd runs

    paths:
        a list containing
        [interactions_dir, ncsd_path, working_dir, output_plotter_path]

    machine:
        name of machine being used for this calculation, e.g. "cedar"

    run:
        boolean, whether or not to run the batch files at the end
    """
    # check manual input
    manual_input_check(man_params, machine, paths)

    # get default parameters
    defaults = Defaults()

    # list of dicts which contain parameters for each run
    list_of_dicts = prepare_input(man_params)
    # creates directories with runnable batch files
    batch_paths = create_dirs(defaults, list_of_dicts, paths, machine)

    # run all batch paths if wanted
    if run:
        print("running all batch files")
        for batch_path in batch_paths:
            if machine == "local":
                system(". "+batch_path)
            elif machine == "cedar":
                system("sbatch "+batch_path)
            elif machine == "summit":
                system("bsub "+batch_path)
            else:
                raise ValueError("Invalid machine!")
    print("done!")
