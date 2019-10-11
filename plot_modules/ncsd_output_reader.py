#!/usr/bin/env python3
"""
-contains function which gets important data from ncsd output files

- also contains a couple small functions for getting names of nuclei
"""

from os.path import split


def element_name(Z):
    element_names = {
        1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O",
        9: "F", 10: "Ne", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",
        16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 21: "Sc", 22: "Ti",
        23: "V", 24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu",
        30: "Zn"
    }   # that's probably all we need, right? Hope so anyway
    return element_names[Z]


def nucleus_name(Z, N):
    """returns the name of a nucleus in 'Li8' style"""
    return element_name(Z) + str(Z+N)


def read_ncsd_output(filename):
    """
    eventually returns a dict of the form:

    dict = {
        "element_name": "B",
        "nucleus_name": "B11",
        "Z": 5,
        "N": 6,
        "Z_plus_N": 11,
        "n_states": 10,
        "filename": "filename",
        "interaction_name": "n3lo-NN3Nlnl-srg2.0"
        "calculated_spectrum": {
            Nmax: {
                state_num: [
                    J(angular momentum),
                    Repetition,
                    Parity,
                    Energy
                ],
            }
        }
    }
    """
    print("reading data from NCSD output file "+filename)
    # read file data
    with open(filename, "r+") as ncsd_output:
        lines = ncsd_output.readlines()

    data_dict = {
        "calculated_spectrum": {},
        "expt_spectrum": {},
        "filename": filename
    }

    # get interaction name from filename
    int_name = split(filename)[1].split("_")[1]
    data_dict["interaction_name"] = int_name

    repetitions = {}

    # get important bits from file
    for line in lines:
        words = line.split()
        if "Z" in line and "N" in line and "hbar" in line:
            Z = int(words[2])
            N = int(words[5])
            data_dict["Z"] = Z
            data_dict["N"] = N
            data_dict["Z_plus_N"] = Z + N
            data_dict["element_name"] = element_name(Z)
            data_dict["nucleus_name"] = nucleus_name(Z, N)
            # the [1:] cuts off an equal sign at the start
            data_dict["hbar_omega"] = float(words[7][1:])
        if "Number of shells" in line:
            n_states = int(words[4])
            data_dict["n_states"] = n_states
        if "Nmax=  " in line:
            # grab Nmax
            words = line.split()
            Nmax = int(words[3])
            data_dict["calculated_spectrum"][Nmax] = {}
            # calculate parity
            # "natural"
            if (N+Z)%2 == 1 and Nmax%2 == 0:
                parity = 1  # negative parity
            elif (N+Z)%2 == 0 and Nmax%2 == 0:
                parity = 0  # positive
            # "unnatural"
            elif (N+Z)%2 == 1 and Nmax%2 == 1:
                parity = 0
            elif (N+Z)%2 == 0 and Nmax%2 == 1:
                parity = 1
            else:
                raise ValueError("What have you done?")
            # empty repetitions dict
            repetitions = {}
        if "State #" in line:
            # parse line for J, repetition, parity, E
            line = line[len("State #")+1:] # to make word separation easier            
            try:
                state_num = int(words[2])
                energy = float(words[5])
                angular_momentum = round(float(words[8]) * 2)  # x2 so it's an int
                isospin = round(float(words[11]) * 2)  # x2 so it's an int                
            except ValueError:
                # for states >= 10, it's written as #10 not # 9
                # I assume we never need to deal with states >=100
                state_num = int(words[1][1:])
                energy = float(words[4])
                angular_momentum = round(float(words[7]) * 2)  # x2 so it's an int
                isospin = round(float(words[10]) * 2)  # x2 so it's an int
            # write down energy relative to state zero energy
            if state_num == 1: # I assume the first state must have already come
                state_1_energy = energy
            # make energy relative
            energy = energy - state_1_energy
            # get repetition number
            if angular_momentum in repetitions.keys():
                repetitions[angular_momentum] += 1
            else:
                repetitions[angular_momentum] = 1
            repetition = repetitions[angular_momentum]
            data_dict["calculated_spectrum"][Nmax][state_num] = [
                angular_momentum, repetition, parity, energy]
    return data_dict


def read_all_ncsd_output(real_paths):
    """runs read_ncsd_output for many paths, tries to merge data"""
    # get data from first file
    overall_data = read_ncsd_output(real_paths[0])

    for path in real_paths[1:]:
        new_data = read_ncsd_output(path)
        # check if okay to merge
        for attr in ["Z", "N", "n_states", "interaction_name"]:
            if new_data[attr] != overall_data[attr]:
                raise ValueError("Data from two files don't match!")
        # merge
        spec = "calculated_spectrum"
        for key in new_data[spec].keys():
            overall_data[spec][key] = new_data[spec][key]
    return overall_data
