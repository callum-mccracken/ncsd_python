"""
A module for calculating parameters to be used in ncsd run files,
using manual parameters.
"""

import os
from . import data_structures
from .formats_multi import kappa_rename_format, potential_end_bit_format


def Nmin_HO(Z):
    """Helper function for Ngs_func"""
    N = 0
    Zrem = Z - 2
    Nmin = 0
    while Zrem > 0:
        N = N + 1
        Nmin = Nmin + N * min((N+1) * (N+2), Zrem)
        Zrem = Zrem - (N+1) * (N+2)
    return Nmin


def Ngs_func(Z, N):
    """
    Calculates Ngs: number of excitations in the ground state

    I didn't want to call this function Ngs, so I could use
    Ngs as a variable later.
    """
    return Nmin_HO(Z) + Nmin_HO(N)


def nucleus(Z, N):
    """Returns the name of a nucleus in 'Li8' style"""
    element_name = {
        1: "H",
        2: "He",
        3: "Li",
        4: "Be",
        5: "B",
        6: "C",
        7: "N",
        8: "O",
        9: "F",
        10: "Ne",
        11: "Na",
        12: "Mg",
        13: "Al",
        14: "Si",
        15: "P",
        16: "S",
        17: "Cl",
        18: "Ar",
        19: "K",
        20: "Ca",
        21: "Sc",
        22: "Ti",
        23: "V",
        24: "Cr",
        25: "Mn",
        26: "Fe",
        27: "Co",
        28: "Ni",
        29: "Cu",
        30: "Zn"
    }   # that's probably all we need, right? Hope so anyway

    return element_name[Z] + str(Z+N)


def calc_params(run_dir, paths, man_params, default_params, machine):
    """
    Calculates all parameters needed for an NCSD run, using manual
    parameters input by the user (man_params) and a set of default_params.

    - run_dir = the directory where the files will be stored
    - paths = list containing
      - interaction files directory
      - path to ncsd executable
      - path to output_plotter module

    Returns MFDPParams, BatchParams objects,
    to be written into a folder for an NCSD run.
    """
    int_dir = paths[0]
    ncsd_path = paths[1]
    output_plotter_filepath = paths[3]
    # for convenience of typing let's make a couple smaller variable names:
    m = man_params
    d = default_params

    Z, N = m.ZN

    # first get nucleus name
    nucleus_name = nucleus(Z, N)

    # we'll need Ngs later:
    Ngs = Ngs_func(Z, N)
    # might as well calculate Nhw now too
    Nhw = m.Nmax_max + Ngs
    nhw_min = m.Nmax_IT + Ngs

    # some harder formatting things to calculate now:
    # output file convention: nucleus_potential_Nmax0-8.freq_IT
    output_file = nucleus_name + "_" + m.potential_name + "_" \
        + "Nmax" + str(m.Nmax_min) + "-" + str(m.Nmax_max) \
        + "." + str(m.hbar_omega)
    if nhw_min <= Nhw:
        output_file += "_IT"

    # string that sets limits on how many nuclei can occupy certain shells
    occupation_string = ""
    for Nshell in range(m.N_1max + 1):
        if Nshell == 0:
            p = 2
            n = 2
            p_plus_n = 4
        elif Nshell == 1:
            p = 6
            n = 6
            p_plus_n = 12
        else:
            # N = shell index starting from 0, m.N = number of neutrons
            p = Z if Z * N <= Nhw else int(Nhw / Nshell)
            n = N if N * N <= Nhw else int(Nhw / Nshell)
            p_plus_n = (Z + N) if (Z + N) * Nshell <= Nhw else int(Nhw / Nshell)
        p = min(p, (N+1)*(N+2))
        n = min(n, (N+1)*(N+2))
        p_plus_n = min(p_plus_n, 2*(N+1)*(N+2))
        line = f" 0 {p}  0 {n}  0 {p_plus_n}  ! N={Nshell}"
        if Nshell != m.N_1max:
            line += "\n"
        occupation_string += line

    # make paths for interaction filess, we'll make these relative paths later
    two_path = os.path.join(int_dir, m.two_body_interaction)
    three_path = os.path.join(int_dir, m.three_body_interaction)
    # remove _comp at the end of 3-body in case it was given by accident
    if three_path[-5:] == "_comp":
        three_path = three_path[:-5]

    # now put all these parameters in a convenient container
    mfdp_parameters = data_structures.MFDPParams(
        # can calculate easily from min_params:
        two_body_interaction=two_path,
        two_body_file_type=m.two_body_interaction[5],
        Z=Z,
        N=N,
        hbar_omega=m.hbar_omega,
        Nhw=Nhw,
        N_1max=m.N_1max,
        N_12max=m.N_12max,
        parity=0 if (Nhw % 2 == 0) else 1,
        total_2Jz=0 if ((Z + N) % 2 == 0) else 1,
        interaction_type=m.interaction_type,
        n_states=m.n_states,
        iterations_required=m.iterations_required,
        irest=m.irest,
        nhw0=m.Nmax_min + Ngs,
        nhw_min=nhw_min,
        nhw_restart=m.nhw_restart,
        kappa_points=m.kappa_points,
        kappa_vals=m.kappa_vals,
        kappa_restart=m.kappa_restart,
        three_body_interaction=three_path,
        N_123max=m.N_123max,
        saved_pivot=m.saved_pivot,
        rmemavail=m.mem,
        # copied from read_params
        N_min=d.N_min,
        iham=d.iham,
        iclmb=d.iclmb,
        strcm=d.strcm,
        major=d.major,
        nshll=d.nshll,
        nsets=d.nsets,
        min_nesp=d.min_nesp,
        nskip=d.nskip,
        iset1=d.iset1,
        ki=d.ki,
        kf=d.kf,
        gs_energy=d.gs_energy,
        igt=d.igt,
        nhme=d.nhme,
        cmin=d.cmin,
        convergence_delta=d.convergence_delta,
        eff_charge_p=d.eff_charge_p,
        eff_charge_n=d.eff_charge_n,
        glp=d.glp,
        gln=d.gln,
        gsp=d.gsp,
        gsn=d.gsn,
        # more complex calculations done earlier
        occupation_string=occupation_string,
        output_file=output_file
    )

    # Now do the batch file's bottom section to do with renaming files.
    non_IT_Nmax = ""
    for i in range(m.Nmax_min, m.Nmax_IT, 2):
        if i <= m.Nmax_max:
            non_IT_Nmax += str(i)+" "
    IT_Nmax = ""
    for i in range(m.Nmax_IT, m.Nmax_max + 1, 2):
        if i <= m.Nmax_max:
            IT_Nmax += str(i)+" "

    # a few functions for converting kappa values to the formats we want
    def kappa_D(kappa_given):
        """Reformats a kappa value, e.g. 2.0 --> 0.200D-04"""
        kappa_e4 = kappa_given * pow(10, -4)
        kappa_scientific = "%.2E" % (kappa_e4)
        kappa_D = kappa_scientific.replace("E", "D")
        [front, back] = kappa_D.split("D")
        converted_front = '0.' + str(int(float(front) * 100))
        kappa = converted_front + "D" + back
        return kappa

    def kappa_em(kappa_given):
        """Reformats a kappa value, e.g. 2.0 --> 2em5"""
        return str(int(kappa_given)) + "em" + "5"

    # get numerical kappa values and create mv lines from those
    kappa_vals = map(float, m.kappa_vals.split())

    kappa_rename = ""
    for i, kappa in enumerate(kappa_vals):
        if i >= m.kappa_points:  # just in case we have more values than needed
            break
        # add each mv line
        kappa_rename += kappa_rename_format.format(
            kappa_D=kappa_D(kappa), kappa_em=kappa_em(kappa)) + "\n"

    potential_end = ""
    if IT_Nmax != "":
        potential_end = potential_end_bit_format.format(
            IT_Nmax=IT_Nmax, kappa_rename=kappa_rename)

    # calculate time
    days, hours, minutes = map(int, m.time.split())
    cedar_time = "{}-{:02}:{:02}".format(days, hours, minutes)
    summit_time = "{}:{:02}".format(days*24 + hours, minutes)

    # organize parameters for batch scripts
    if machine == "local":
        batch_parameters = data_structures.LocalBatchParams(
            run_directory=run_dir,
            potential=m.potential_name,
            nucleus_name=nucleus_name,
            hbar_omega=int(m.hbar_omega),
            suffix="_"+str(m.n_states)+"st",
            Ngs=Ngs,
            ncsd_path=ncsd_path,
            non_IT_Nmax=non_IT_Nmax,
            potential_end_bit=potential_end,
            output_file=output_file,
            output_plotter=output_plotter_filepath
        )
    elif machine == "cedar":
        batch_parameters = data_structures.CedarBatchParams(
            run_directory=run_dir,
            account="rrg-navratil",
            nodes=m.n_nodes,
            tasks_per_node=d.tasks_per_node,
            mem_per_core=int(m.mem),
            time=cedar_time,
            output="ncsd-%J.out",
            potential=m.potential_name,
            nucleus_name=nucleus_name,
            hbar_omega=int(m.hbar_omega),
            suffix="_"+str(m.n_states)+"st",
            Ngs=Ngs,
            ncsd_path=ncsd_path,
            non_IT_Nmax=non_IT_Nmax,
            potential_end_bit=potential_end,
            output_file=output_file,
            output_plotter=output_plotter_filepath
        )
    elif machine == "summit":
        batch_parameters = data_structures.SummitBatchParams(
            run_directory=run_dir,
            account="nph123",
            nnodes=m.n_nodes,
            time=summit_time,
            resource_sets=6 * m.n_nodes,
            output="ncsd-run_"+nucleus_name+".out",
            potential=m.potential_name,
            nucleus_name=nucleus_name,
            hbar_omega=int(m.hbar_omega),
            suffix="_"+str(m.n_states)+"st",
            Ngs=Ngs,
            ncsd_path=ncsd_path,
            non_IT_Nmax=non_IT_Nmax,
            potential_end_bit=potential_end,
            output_file=output_file,
            output_plotter=output_plotter_filepath
        )
    else:
        raise ValueError("Invalid machine: "+machine)
    return mfdp_parameters, batch_parameters
