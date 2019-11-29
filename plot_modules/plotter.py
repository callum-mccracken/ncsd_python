"""Contains functions to help plot (or otherwise export) data."""
from . import formats

import os
import numpy as np
import matplotlib.pyplot as plt

"""
data must be of the form:
    data = {
        "skip_Nmax" = [],
        "max_state" = 10,
        "nucleus_name": "B11",
        "Z": 5,
        "N": 6,
        "n_states": 10,    
        "element_name": "Li",
        "Z_plus_N": "8",
        "interaction_name": "some_name",
        "filename": "filename,
        "calculated_spectrum": {
            <Nmax>: {
                <state_num>: [
                    angular_momentum: float,
                    repetition: integer,
                    parity: integer,
                    energy: float
                ]
            }
            ...
        }
        "expt_spectrum": {
            "Expt": {
                <state_num>: [
                    angular_momentum: float,
                    repetition: integer,
                    parity: integer,
                    energy: float
                ]
            }
        }
    }
"""

def write_xmgrace(input_data, save_dir, grace_plotter_path):
    """creates a file which can be used by xmgrace"""
    # NOTE: "ex" prefix --> excitation energies (as opposed to binding energies)

    # let's create the datasets and axis labels first
    c_spectrum = input_data["calculated_spectrum"]
    e_spectrum = input_data["expt_spectrum"]
    
    # variables to contain numerical data
    data_string = ""
    ex_data_string = ""

    # string containing axis labels
    axis_labels = ""

    # the maximum state number that will be plotted
    max_state = input_data["max_state"]
    
    for Nmax in sorted(c_spectrum.keys()):
        if Nmax in input_data["skip_Nmax"]:
            continue
        title = formats.xmgrace_Nmax_title_format.format(Nmax=Nmax)
        # lines = string to hold the lines containing
        # [j, repetition, parity, energy]
        # that we're going to generate
        lines = ""
        ex_lines = ""
        
        # this stores the ground state energy, for calculating relative energies
        state_1_energy = c_spectrum[Nmax][1][3]
        for state_num in sorted(c_spectrum[Nmax].keys()):
            if state_num > max_state:
                continue
            j, repetition, parity, energy = c_spectrum[Nmax][state_num]
            ex_energy = energy - state_1_energy
            lines += formats.xmgrace_data_line_format.format(
                j=j,repetition=repetition, parity=parity, energy=energy)
            ex_lines += formats.xmgrace_data_line_format.format(
                j=j, repetition=repetition, parity=parity, energy=ex_energy)

        data_string += formats.xmgrace_dataset_format.format(
            title=title, lines=lines)
        ex_data_string += formats.xmgrace_dataset_format.format(
            title=title, lines=ex_lines)
        axis_labels += formats.xmgrace_axis_label_line.format(Nmax=Nmax)
    
    # experimental data, only 1 dataset.
    title = "Expt"
    lines = ""
    ex_lines = ""
    for state_num in sorted(e_spectrum[title].keys()):
        if state_num > max_state:
            continue
        j, repetition, parity, energy = c_spectrum[Nmax][state_num]
        # note: state_1_energy is the ground state energy of the largest Nmax
        # (since we haven't changed it since the end of the last loop)
        ex_energy = energy - state_1_energy
        lines += formats.xmgrace_data_line_format.format(
            j=j, repetition=repetition, parity=parity, energy=energy)
        ex_lines += formats.xmgrace_data_line_format.format(
            j=j, repetition=repetition, parity=parity, energy=ex_energy)
    data_string += formats.xmgrace_dataset_format.format(
        title=title, lines=lines)
    ex_data_string += formats.xmgrace_dataset_format.format(
        title=title, lines=ex_lines)
    axis_labels += title
    
    # now a couple final things
    num_spectra = len(c_spectrum.keys()) + len(e_spectrum.keys())
    num_spectra = num_spectra - len(input_data["skip_Nmax"])
    num_states = len(e_spectrum["Expt"].keys())
    if input_data["max_state"] < num_states:
        num_states = input_data["max_state"] 
    Z_plus_N = input_data["Z_plus_N"]
    element = input_data["element_name"]
    interaction =input_data["interaction_name"]
    
    # then save absolute energies and relative ones
    filename = os.path.split(input_data["filename"])[-1]
    filename = filename[:filename.index("_Nmax")]+'_spectra_vs_Nmax.grdt'
    ex_filename = filename[:filename.index(".grdt")]+'_excited.grdt'
    save_path = os.path.join(save_dir, filename)
    ex_save_path = os.path.join(save_dir, ex_filename)
    # write binding energy file
    with open(save_path, "w+") as open_file:
        open_file.write(
            formats.xmgrace_format.format(
                num_spectra_plus_2 = num_spectra + 2,
                num_states = num_states,
                num_spectra = num_spectra,
                num_plots = 1,
                Z_plus_N = Z_plus_N,
                element = element,
                interaction_name = interaction,
                axis_labels = axis_labels,
                data = data_string
            )
        )
    # write excited file
    with open(ex_save_path, "w+") as open_file:
        open_file.write(
            formats.xmgrace_format.format(
                num_spectra_plus_2 = num_spectra + 2,
                num_states = num_states,
                num_spectra = num_spectra,
                num_plots = 1,
                Z_plus_N = Z_plus_N,
                element = element,
                interaction_name = interaction,
                axis_labels = axis_labels,
                data = ex_data_string
            )
        )
    
    # now call the grace_spectra_plotter.exe file
    os.system(grace_plotter_path + " " + save_path)
    os.system(grace_plotter_path + " " + save_path + " -excited")
    plot = False
    if plot:
        # and actually use xmgrace to plot
        agr_path = save_path[:save_path.index(".grdt")] + ".agr"
        os.system("xmgrace " + agr_path)
        ex_agr_path = save_path[:ex_save_path.index(".grdt")] + ".agr"
        os.system("xmgrace " + ex_agr_path)


def write_csv(input_data, save_dir):
    """creates a csv file containing useful data, for easier parsing later"""
    # write titles
    file_string = ",".join(
        ["Title", "StateNum", "J", "repetition", "Parity", "Energy"]) + "\n"
    
    # let's create the datasets and axis labels first
    c_spectrum = input_data["calculated_spectrum"]
    e_spectrum = input_data["expt_spectrum"]
    
    # calculated data
    max_state = input_data["max_state"]
    for Nmax in sorted(c_spectrum.keys()):
        if Nmax in input_data["skip_Nmax"]:
            continue
        title = "Nmax"+str(Nmax)
        lines = ""
        for state_num in sorted(c_spectrum[Nmax].keys()):
            if state_num > max_state:
                continue            
            j=str(c_spectrum[Nmax][state_num][0])
            repetition=str(c_spectrum[Nmax][state_num][1])
            parity=str(c_spectrum[Nmax][state_num][2])
            energy=str(c_spectrum[Nmax][state_num][3])
            lines +=  ",".join(
                [title, str(state_num), j, repetition, parity, energy]) + "\n"
        file_string += lines

    # experimental data, only 1 dataset
    title = "Expt"
    lines = ""
    for state_num in sorted(e_spectrum[title].keys()):
        if state_num > max_state:
            continue
        j=str(e_spectrum[title][state_num][0])
        repetition=str(e_spectrum[title][state_num][1])
        parity=str(e_spectrum[title][state_num][2])
        energy=str(e_spectrum[title][state_num][3])
        lines += ",".join(
            [title, str(state_num), j, repetition, parity, energy]) + "\n"
    file_string += lines
    # now save
    filename = os.path.split(input_data["filename"])[-1]
    filename = filename[:filename.index("_Nmax")]+'_spectra_vs_Nmax.csv'
    with open(os.path.join(save_dir, filename), "w+") as open_file:
        open_file.write(file_string)


def matplotlib_plot(input_data, save_dir):
    """makes a matplotlib style plot of our data"""

    energy_datasets = []
    axis_labels = []
    
    # let's create the datasets and axis labels first
    c_spectrum = input_data["calculated_spectrum"]
    e_spectrum = input_data["expt_spectrum"]
    
    # calculated data
    max_state = input_data["max_state"]
    for Nmax in sorted(c_spectrum.keys()):
        if Nmax in input_data["skip_Nmax"]:
            continue
        axis_labels.append(str(Nmax)+"$\\hbar \\omega$")
        dataset = []
        for state_num in sorted(c_spectrum[Nmax].keys()):
            if state_num > max_state:
                continue            
            energy=float(c_spectrum[Nmax][state_num][3])
            dataset.append(energy)
        energy_datasets.append(dataset)

    # experimental data, only 1 dataset
    axis_labels.append("Expt")
    dataset = []
    for state_num in sorted(e_spectrum["Expt"].keys()):
        if state_num > max_state:
            continue
        energy=float(e_spectrum["Expt"][state_num][3])
        dataset.append(energy)
    energy_datasets.append(dataset)

    energy_datasets = np.array(energy_datasets)
    plot_arr = energy_datasets.transpose()

    # create figure
    f = plt.figure()
    ax = f.add_subplot(111)
    plt.ylabel("$E_x$ [MeV]")
    
    # play around with ticks
    #ax.yaxis.tick_right()
    ax.yaxis.set_tick_params(right=True, direction="in")
    ax.xaxis.set_tick_params(length=0)
    padding = 2
    e_min = int(round(np.amin(plot_arr)) - padding)
    e_max = int(round(np.amax(plot_arr)) + padding)
    plt.yticks(range(e_min, e_max+1))
    plt.xticks(np.arange(0.5, 2*len(plot_arr[0])-0.5, 2.0))
    ax.set_xticklabels(axis_labels)
    
    # TODO: do something with the title / element name and line labels
    # TODO: make this prettier in general

    for line in plot_arr:
        # pick a colour for the line
        colour = np.random.rand(3,)
        for i, value in enumerate(line):
            # plot in a bit of a strange way, so it keeps the xmgrace format.
            # solid line for data point
            plt.plot([2*i,2*i+1],[value, value],
                c=colour, linestyle="solid")
            if i != len(line) - 1:
                # dotted line to show how it changes
                next_val = line[i+1]
                plt.plot([2*i+1, 2*i+2], [value, next_val],
                    c=colour, linestyle="dotted")
            # othewise we've reached the end of the list    

    # save
    filename = os.path.split(input_data["filename"])[-1]
    filename = filename[:filename.index("_Nmax")]+'_spectra_vs_Nmax.png'
    plt.savefig(os.path.join(save_dir, filename))


def export_data(data, save_dir, grace_plotter_path, out_type="xmgrace"):
    # gives a bunch of different ways to produce output
    if out_type == "xmgrace":
        write_xmgrace(data, save_dir, grace_plotter_path)
    elif out_type == "csv":
        write_csv(data, save_dir)
    elif out_type == "matplotlib":
        matplotlib_plot(data, save_dir)
    else:
        raise ValueError("The output type "+out_type+" is not supported yet")    