"""
Contains functions to help plot (or otherwise export) data.

When calling most of these functions, data must be of the form::

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
                    repetition: int,
                    parity: int,
                    energy: float
                ]
            }
            ...
        }
        "expt_spectrum": {
            "Expt": {
                <state_num>: [
                    angular_momentum: float,
                    repetition: int,
                    parity: int,
                    energy: float
                ]
            }
        }
    }

"""
from . import formats_plot

import os
import numpy as np
import matplotlib.pyplot as plt


def write_xmgrace(input_data, save_dir, grace_plotter_path):
    """
    Creates a .agr file which can be plotted by xmgrace,
    or rather just a .grdt file,
    which you have to convert to .agr using grace_plotter.exe.

    input_data:
        a dictionary of the form described at the top of this file

    save_dir:
        path to the directory where we should save the plot

    grace_plotter_path:
        the path to the xmgrace plotter executable file
    """
    # NOTE: "ex" prefix --> excitation energies (not binding energies)

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
        title = formats_plot.xmgrace_Nmax_title_format.format(Nmax=Nmax)
        # lines = string to hold the lines containing
        # [j, repetition, parity, energy]
        # that we're going to generate
        lines = ""
        ex_lines = ""

        # stores the ground state energy, for calculating relative energies
        state_1_energy = c_spectrum[Nmax][1][3]
        for state_num in sorted(c_spectrum[Nmax].keys()):
            if state_num > max_state:
                continue
            j, repetition, parity, energy = c_spectrum[Nmax][state_num]
            ex_energy = energy - state_1_energy
            lines += formats_plot.xmgrace_data_line_format.format(
                j=j, repetition=repetition, parity=parity, energy=energy)
            ex_lines += formats_plot.xmgrace_data_line_format.format(
                j=j, repetition=repetition, parity=parity, energy=ex_energy)

        data_string += formats_plot.xmgrace_dataset_format.format(
            title=title, lines=lines)
        ex_data_string += formats_plot.xmgrace_dataset_format.format(
            title=title, lines=ex_lines)
        axis_labels += formats_plot.xmgrace_axis_label_line.format(Nmax=Nmax)

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
        lines += formats_plot.xmgrace_data_line_format.format(
            j=j, repetition=repetition, parity=parity, energy=energy)
        ex_lines += formats_plot.xmgrace_data_line_format.format(
            j=j, repetition=repetition, parity=parity, energy=ex_energy)
    data_string += formats_plot.xmgrace_dataset_format.format(
        title=title, lines=lines)
    ex_data_string += formats_plot.xmgrace_dataset_format.format(
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
    interaction = input_data["interaction_name"]

    # then save absolute energies and relative ones
    filename = os.path.split(input_data["filename"])[-1]
    filename = filename[:filename.index("_Nmax")]+'_spectra_vs_Nmax.grdt'
    ex_filename = filename[:filename.index(".grdt")]+'_excited.grdt'
    save_path = os.path.join(save_dir, filename)
    ex_save_path = os.path.join(save_dir, ex_filename)
    # write binding energy file
    with open(save_path, "w+") as open_file:
        open_file.write(
            formats_plot.xmgrace_format.format(
                num_spectra_plus_2=num_spectra + 2,
                num_states=num_states,
                num_spectra=num_spectra,
                num_plots=1,
                Z_plus_N=Z_plus_N,
                element=element,
                interaction_name=interaction,
                axis_labels=axis_labels,
                data=data_string
            )
        )
    # write excited file
    with open(ex_save_path, "w+") as open_file:
        open_file.write(
            formats_plot.xmgrace_format.format(
                num_spectra_plus_2=num_spectra + 2,
                num_states=num_states,
                num_spectra=num_spectra,
                num_plots=1,
                Z_plus_N=Z_plus_N,
                element=element,
                interaction_name=interaction,
                axis_labels=axis_labels,
                data=ex_data_string
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
    """
    Creates a csv file containing useful data, for easier parsing later

    input_data:
        a dictionary of the form described at the top of this file

    save_dir:
        path to a diretory where we should save the csv file
    """
    # write titles
    file_string = ",".join(
        ["Title", "StateNum", "J", "repetition", "Parity", "AbsEnergy"]) + "\n"
    ex_file_string = ",".join(
        ["Title", "StateNum", "J", "repetition", "Parity", "ExcEnergy"]) + "\n"

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
        ex_lines = ""
        state_1_energy = c_spectrum[Nmax][1][3]
        for state_num in sorted(c_spectrum[Nmax].keys()):
            if state_num > max_state:
                continue
            j = str(c_spectrum[Nmax][state_num][0])
            repetition = str(c_spectrum[Nmax][state_num][1])
            parity = str(c_spectrum[Nmax][state_num][2])
            energy = str(c_spectrum[Nmax][state_num][3])
            ex_energy = energy - state_1_energy
            lines += ",".join(
                [title, str(state_num), j, repetition, parity, energy]) + "\n"
            ex_lines += ",".join(
                [title, str(state_num), j, repetition, parity, ex_energy]) + "\n"

        file_string += lines
        ex_file_string += ex_lines

    # experimental data, only 1 dataset
    title = "Expt"
    lines = ""
    ex_lines = ""
    for state_num in sorted(e_spectrum[title].keys()):
        if state_num > max_state:
            continue
        j = str(e_spectrum[title][state_num][0])
        repetition = str(e_spectrum[title][state_num][1])
        parity = str(e_spectrum[title][state_num][2])
        energy = str(e_spectrum[title][state_num][3])
        ex_energy = energy - state_1_energy
        lines += ",".join(
            [title, str(state_num), j, repetition, parity, energy]) + "\n"
        ex_lines += ",".join(
            [title, str(state_num), j, repetition, parity, ex_energy]) + "\n"
    file_string += lines
    ex_file_string += ex_lines
    # now save
    filename = os.path.split(input_data["filename"])[-1]
    filename = filename[:filename.index("_Nmax")]+'_spectra_vs_Nmax.csv'
    with open(os.path.join(save_dir, filename), "w+") as open_file:
        open_file.write(file_string)
    ex_filename = filename.replace(".csv", "_excited.csv")
    with open(os.path.join(save_dir, ex_filename), "w+") as open_file:
        open_file.write(ex_file_string)


def matplotlib_plot(input_data, save_dir):
    """
    Makes a matplotlib style plot of our data, saves as .png, .svg

    input_data:
        a dictionary of the form described at the top of this file

    save_dir:
        path to a diretory where we should save our plots
    """

    energies = []
    ex_energies = []
    axis_labels = []
    line_labels = []

    # let's create the datasets and axis labels first
    c_spectrum = input_data["calculated_spectrum"]
    e_spectrum = input_data["expt_spectrum"]

    # calculated data
    max_state = input_data["max_state"]
    for Nmax in sorted(c_spectrum.keys()):
        if Nmax in input_data["skip_Nmax"]:
            continue
        axis_labels.append(str(Nmax)+"$\\hbar \\omega$")
        e_list = []
        ex_e_list = []
        state_1_energy = c_spectrum[Nmax][1][3]
        for state_num in sorted(c_spectrum[Nmax].keys()):
            if state_num > max_state:
                continue
            angular_momentum = int(c_spectrum[Nmax][state_num][0])
            repetition = int(c_spectrum[Nmax][state_num][1])
            parity = int(c_spectrum[Nmax][state_num][2])
            energy = float(c_spectrum[Nmax][state_num][3])
            ex_energy = energy - state_1_energy
            title = f"${angular_momentum}^{parity}$"
            e_list.append(energy)
            ex_e_list.append(ex_energy)
        energies.append(e_list)
        ex_energies.append(ex_e_list)

    # experimental data, only 1 dataset
    axis_labels.append("Expt")
    e_list = []
    ex_e_list = []
    for state_num in sorted(e_spectrum["Expt"].keys()):
        if state_num > max_state:
            continue
        angular_momentum = int(e_spectrum["Expt"][state_num][0])
        parity = e_spectrum["Expt"][state_num][2]
        parity = "+" if parity == "0" else "-"
        energy = float(e_spectrum["Expt"][state_num][3])
        ex_energy = energy - state_1_energy
        title = f"${angular_momentum}^{parity}$"
        line_labels.append(title)
        e_list.append(energy)
        ex_e_list.append(ex_energy)
    energies.append(e_list)
    ex_energies.append(ex_e_list)

    energies = np.array(energies)
    ex_energies = np.array(ex_energies)
    plot_arr = energies.transpose()
    ex_plot_arr = energies.transpose()

    # create figure
    f = plt.figure()
    ax = f.add_subplot(111)
    plt.ylabel("$E$ [MeV]")

    # play around with ticks
    # ax.yaxis.tick_right()
    ax.yaxis.set_tick_params(right=True, direction="in")
    ax.xaxis.set_tick_params(length=0)
    plt.xticks(np.arange(0.5, 2*len(plot_arr[0])-0.5, 2.0))
    ax.set_xticklabels(axis_labels)
    plt.xlim(-1, 2*len(plot_arr[0])+1)

    for num, line in enumerate(plot_arr):
        # pick a colour for the line
        colour = np.random.rand(3,)
        for i, value in enumerate(line):
            # plot in a bit of a strange way, so it keeps the xmgrace format.
            # solid line for data point
            plt.plot([2*i, 2*i+1], [value, value],
                     c=colour, linestyle="solid")
            if i != len(line) - 1:
                # dotted line to show how it changes
                next_val = line[i+1]
                plt.plot([2*i+1, 2*i+2], [value, next_val],
                         c=colour, linestyle="dotted")
            # othewise we've reached the end of the list
        plt.text(2*len(line) - 0.5, line[-1], line_labels[num])

    # save the plot
    filename = os.path.split(input_data["filename"])[-1]
    plt_title = filename.split("_")[0] + " Bound States"
    plt.title(plt_title)
    filename = filename[:filename.index("_Nmax")]+'_spectra_vs_Nmax'
    plt.savefig(os.path.join(save_dir, filename+".png"))
    plt.savefig(os.path.join(save_dir, filename+".svg"))

    # same code but for excitation energies

    plt.cla(); plt.clf()
    # create figure
    f = plt.figure()
    ax = f.add_subplot(111)
    plt.ylabel("$E_x$ [MeV]")

    # play around with ticks
    # ax.yaxis.tick_right()
    ax.yaxis.set_tick_params(right=True, direction="in")
    ax.xaxis.set_tick_params(length=0)
    plt.xticks(np.arange(0.5, 2*len(ex_plot_arr[0])-0.5, 2.0))
    ax.set_xticklabels(axis_labels)
    plt.xlim(-1, 2*len(ex_plot_arr[0])+1)

    for num, line in enumerate(ex_plot_arr):
        # pick a colour for the line
        colour = np.random.rand(3,)
        for i, value in enumerate(line):
            # plot in a bit of a strange way, so it keeps the xmgrace format.
            # solid line for data point
            plt.plot([2*i, 2*i+1], [value, value],
                     c=colour, linestyle="solid")
            if i != len(line) - 1:
                # dotted line to show how it changes
                next_val = line[i+1]
                plt.plot([2*i+1, 2*i+2], [value, next_val],
                         c=colour, linestyle="dotted")
            # othewise we've reached the end of the list
        plt.text(2*len(line) - 0.5, line[-1], line_labels[num])

    # save the plot
    filename = filename+"_excited"
    plt.title(plt_title)
    plt.savefig(os.path.join(save_dir, filename+".png"))
    plt.savefig(os.path.join(save_dir, filename+".svg"))


def export_data(data, save_dir, grace_plotter_path, out_type="xmgrace"):
    """
    Produces output in a bunch of different ways.

    data:
        a dictionary of the form described at the top of this file

    save_dir:
        path to a diretory where we should save our plots

    grace_plotter_path:
        the path to the xmgrace plotter executable file

    out_type:
        string, either "xmgrace", "csv", or "matplotlib"
    """
    if out_type == "xmgrace":
        write_xmgrace(data, save_dir, grace_plotter_path)
    elif out_type == "csv":
        write_csv(data, save_dir)
    elif out_type == "matplotlib":
        matplotlib_plot(data, save_dir)
    else:
        raise ValueError("The output type "+out_type+" is not supported yet")
