"""
Contains functions to check various sorts of data.
"""
import sys
from os.path import join, exists
import os
import traceback

def get_int_dir():
    """
    Check if INT_DIR environment variable exists, create it if not.

    Not used by default, but it's an option if you want it.
    """
    try:
        int_dir = os.environ["INT_DIR"]
    except KeyError:
        int_dir = input("Enter the (full path) directory "
                        "where your interactions are stored: ")
        os.system(f"echo 'export INT_DIR=\"{int_dir}\"\n' >> ~/.bash_profile")
        os.system(". ~/.bash_profile")
    return int_dir


def get_work_dir():
    """
    Check if WORK_DIR environment variable exists, create it if not.

    Not used by default, but it's an option if you want it.
    """
    try:
        wrk_dir = os.environ["WORK_DIR"]
    except KeyError:
        wrk_dir = input("Enter the full path to your working directory: ")
        os.system(f"echo 'export WORK_DIR=\"{wrk_dir}\"\n' >> ~/.bash_profile")
        os.system(". ~/.bash_profile")
    return wrk_dir


def manual_input_check(manual_params, machine, paths):
    """
    Checks manual input to ensure it is self-consistent.

    manual_params:
        a data_structures.ManParams object, filled with manual input variables

    machine:
        string designating which machine is being used

    paths:
        list containing: [interactions dir, ncsd exe path, working dir]
    """
    print("checking manual input")
    m = manual_params  # so we don't have to type out manual_params everywhere

    int_dir = paths[0]
    ncsd_path = paths[1]
    working_dir = paths[2]
    # do we have a 3-body interaction?
    three_body = (abs(m.interaction_type) == 3)

    # first check if paths exist
    if not exists(int_dir):
        raise IOError(
            "Interactions directory " + int_dir + " does not exist")
    if not exists(working_dir):
        raise IOError(
            "Working directory " + working_dir + " does not exist")
    f2 = join(int_dir, m.two_body_interaction)
    if not exists(f2):
        raise IOError("Two body file "+f2+" does not exist")
    if three_body:
        f3 = join(int_dir, m.three_body_interaction)
        if not exists(f3):
            raise IOError("Three body file "+f3+" does not exist")
    if not exists(ncsd_path):
        raise IOError("NCSD file "+ncsd_path+" does not exist!")

    # check that parameters make sense
    if not (m.N_12max >= m.N_1max):
        raise ValueError("N_12max must be >= N_1max")
    if three_body:
        if not (m.N_123max >= m.N_12max):
            raise ValueError("N_123max must be >= N_12max")

    # check that parameters match with filenames
    try:
        # TBME file
        tbme_filename = m.two_body_interaction
        last_chunk = tbme_filename.split(".")[-1]
        [hbar_omega_verif_0, other_stuff] = last_chunk.split("_")
        hbar_omega_verif_0 = float(hbar_omega_verif_0)
        # see if str(N_1max) + str(N_1max) == other_stuff
        if other_stuff != str(m.N_1max) + str(m.N_12max):
            print("\nYour TMBE file doesn't seem to match your parameters!")
            print("N_1max = "+str(m.N_1max))
            print("N_12max = "+str(m.N_12max))
            print("TBME filename = "+tbme_filename)
            print("relevant section = "+other_stuff)
            yn = ""
            while yn not in ["y", "n"]:
                yn = input("Do you want to continue? (y/n): ")
            if yn == "y":
                pass
            else:
                sys.exit(0)
        # see if hbar_omega matches
        if hbar_omega_verif_0 != m.hbar_omega:
            print("\nYour TMBE file doesn't seem to match your parameters!")
            print("hbar_omega = "+str(m.hbar_omega))
            print("TBME filename = "+tbme_filename)
            print("hbar_omega from the file is", hbar_omega_verif_0)
            yn = ""
            while yn not in ["y", "n"]:
                yn = input("Do you want to continue? (y/n): ")
            if yn == "y":
                pass
            else:
                sys.exit(0)
    except Exception as e:
        print("Minor error caught while parsing TMBE filename.")
        print("Printing traceback as if it had caused a crash:")
        traceback.print_exc()
        print("TBME filename that caused this error:", tbme_filename)
        print("We assume everything's fine, but double-check!\n")

    if three_body:
        try:
            # three-body file
            three_filename = m.three_body_interaction
            [penultimate_chunk, last_chunk] = three_filename.split(".")[-2:]
            # get hbar_omega
            [hbar_omega_verif_1, other_stuff] = last_chunk.split("_")
            hbar_omega_verif_1 = float(hbar_omega_verif_1)
            # get N_#max variables
            n_maxes = penultimate_chunk.split("_")[-1]
            # see if str(N_1max) + str(N_1max) == other_stuff
            if n_maxes != str(m.N_123max) + str(m.N_12max) + str(m.N_1max):
                print(
                    "\nYour 3-body file doesn't seem "
                    "to match your parameters!")
                print("N_1max = "+str(m.N_1max))
                print("N_12max = "+str(m.N_12max))
                print("N_123max = "+str(m.N_123max))
                print("3-body filename = "+three_filename)
                print("relevant section = "+n_maxes)
                yn = ""
                while yn not in ["y", "n"]:
                    yn = input("Do you want to continue? (y/n): ")
                if yn == "y":
                    pass
                else:
                    sys.exit(0)
            # see if hbar_omega matches
            if hbar_omega_verif_1 != m.hbar_omega:
                print(
                    "\nYour 3-body file doesn't seem "
                    "to match your parameters!")
                print("hbar_omega = "+str(m.hbar_omega))
                print("3-body filename = "+three_filename)
                print("hbar_omega from the file is", hbar_omega_verif_1)
                yn = ""
                while yn not in ["y", "n"]:
                    yn = input("Do you want to continue? (y/n): ")
                if yn == "y":
                    pass
                else:
                    sys.exit(0)
        except Exception as e:
            print("Minor error caught while parsing 3-body filename.")
            print("Printing traceback as if it had caused a crash:")
            traceback.print_exc()
            print("3-body filename that caused the error:", three_filename)
            print("We assume everything's fine, but double-check!\n")

    # check there's at least kappa_points kappa values
    kappa_vals = list(map(float, m.kappa_vals.split()))
    if len(kappa_vals) < m.kappa_points:
        raise ValueError(
            "You must have at least kappa_points kappa values!"
            " kappa_points = "+str(m.kappa_points))

    # and if kappa_points and kappa_vals disagree, make sure they know that
    if len(kappa_vals) > m.kappa_points:
        print(
            "Did you mean to enter "+str(len(kappa_vals)) +
            " values for kappa_min, but set kappa_points to " +
            str(m.kappa_points)+"?")
        user_input = ""
        while user_input not in ["Y", "N"]:
            user_input = input("Enter Y to proceed, N to cancel: ")
        if user_input == "N":
            print("Okay, exiting... Try again!")
            sys.exit(0)

    kr_values = [-1, 1, 2, 3, 4]
    if m.kappa_restart not in kr_values:
        raise ValueError(
            "kappa_restart must be one of" + " ".join(map(str, kr_values)))

    if m.saved_pivot not in ["F", "T"]:
        raise ValueError("saved_pivot must be either T or F")

    if (m.irest == 1 or m.kappa_restart != -1 or m.nhw_restart != -1) \
       and m.saved_pivot == "F":
        raise ValueError("why not use the saved pivot if you're restarting?")

    #  if this function runs, the input passes the test


def check_mfdp_read(mfdp_params):
    """
    Checks to see if mfdp data, read from a file, was ok.

    Note that this is not actually used now, since we don't try to parse
    mfdp files anymore, it's just left around from when we did, in case
    you want to use it for something later.

    mfdp_params:
        a data_structures.MFDPParams object
    """
    print("opening mfdp file, checking data")

    # 3 body interaction?
    three_body = (abs(mfdp_params.interaction_type) == 3)

    # there must be a better way than hard-coding all these indices, right?

    # parse TBME file for parameters
    directories = mfdp_params.two_body_interaction.split("/")
    tbme_filename = directories[-1]  # last one's the actual file
    tbme_type = int(tbme_filename[5])
    last_chunk = tbme_filename.split(".")[-1]
    [hbar_omega_verif_0, other_stuff] = last_chunk.split("_")
    hbar_omega_verif_0 = float(hbar_omega_verif_0)
    N_1max_verif = int(other_stuff[0])
    N_12max_verif = int(other_stuff[1:])

    # parse output file name
    sections = mfdp_params.output_file.split("_")
    the_rest = sections[2]
    dot_sections = the_rest.split(".")
    hbar_omega_verif_1 = float(dot_sections[1])

    # parse 3-body
    if three_body:
        # check 3-body file?
        pass

    message = ""  # adjust error message as needed

    # check obvious things
    if mfdp_params.saved_pivot not in ["F", "T"]:
        message = "saved_pivot must be either T or F"
    if mfdp_params.two_body_file_type != tbme_type:
        message = "TBME type does not match type from TMBE filename"
    if mfdp_params.hbar_omega != hbar_omega_verif_0:
        message = "freq does not match freq from TMBE filename"
    if mfdp_params.hbar_omega != hbar_omega_verif_1:
        message = "freq does not match freq from output filename"
    if mfdp_params.N_1max != N_1max_verif:
        message = "N_1max does not match value from TBME filename"
    if mfdp_params.N_12max != N_12max_verif:
        message = "N_12max does not match value from TBME filename"
    if mfdp_params.eff_charge_p != 1.0:
        message = ("effective charge of proton is 1.0, "
                   "not "+str(mfdp_params.eff_charge_p))
    if mfdp_params.eff_charge_n != 0.0:
        message = ("effective charge of neutron is 0.0, "
                   "not "+str(mfdp_params.eff_charge_n))
    if mfdp_params.glp != 1.0:
        message = "glp is always 1.0, not "+str(mfdp_params.glp)
    if mfdp_params.gln != 0.0:
        message = "gln is always 1.0, not "+str(mfdp_params.gln)
    if mfdp_params.gsp != 5.586:
        message = "gsp is always 5.586, not "+str(mfdp_params.gsp)
    if mfdp_params.gsn != -3.826:
        message = "gsn is always -3.826, not "+str(mfdp_params.gsn)

    # mod 2 checks
    Z, N = mfdp_params.ZN
    if ((Z + N) % 2) == 0:
        if mfdp_params.total_2Jz != 0:
            message = ("Z + N is even, so total_2Jz must be 0, "
                       "not "+str(mfdp_params.total_2Jz))
    else:
        if mfdp_params.total_2Jz != 1:
            message = ("Z + N is odd, so total_2Jz must be 1, "
                       "not "+str(mfdp_params.total_2Jz))

    if mfdp_params.parity != (mfdp_params.Nhw % 2):
        message = "we require parity = Nhw mod 2"
    if mfdp_params.parity != (mfdp_params.nhw0 % 2):
        message = "we require parity = nhw0 mod 2"
    if mfdp_params.parity != (mfdp_params.nhw_min % 2):
        message = "we require parity = nhw_min mod 2"

    # raise last error detected
    if message:
        raise ValueError("Bad template MFDP file: "+message)
