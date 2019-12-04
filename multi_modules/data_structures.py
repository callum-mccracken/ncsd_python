"""
Contains the Params class and subclasses for specific file types.
"""

man_keys = [
    "two_body_interaction",
    "Z",
    "N",
    "hbar_omega",
    "N_1max",
    "N_12max",
    "N_123max",
    "three_body_interaction",
    "Nmax_min",
    "Nmax_max",
    "Nmax_IT",
    "interaction_type",
    "n_states",
    "iterations_required",
    "irest",
    "nhw_restart",
    "kappa_points",
    "kappa_vals",
    "kappa_restart",
    "saved_pivot",
    "mem",
    "time",
    "potential_name",
    "n_nodes"]
"""all manual input variables"""

local_batch_keys = [
    "run_directory",
    "potential",
    "nucleus_name",
    "hbar_omega",
    "suffix",
    "Ngs",
    "ncsd_path",
    "non_IT_Nmax",
    "potential_end_bit",
    "output_file",
    "output_plotter"
    ]
"""all variables for batch file on local machine"""

cedar_batch_keys = [
    "run_directory",
    "account",
    "nodes",
    "tasks_per_node",
    "mem_per_core",
    "mem",
    "time",
    "output",
    "potential",
    "nucleus_name",
    "hbar_omega",
    "suffix",
    "Ngs",
    "ncsd_path",
    "non_IT_Nmax",
    "potential_end_bit",
    "output_file",
    "output_plotter"
    ]
"""all variables for batch file on cedar"""

summit_batch_keys = [
    "run_directory",
    "account",
    "nnodes",
    "resource_sets",
    "time",
    "output",
    "potential",
    "nucleus_name",
    "hbar_omega",
    "suffix",
    "Ngs",
    "ncsd_path",
    "non_IT_Nmax",
    "potential_end_bit",
    "output_file",
    "output_plotter"
    ]
"""all variables for batch file on summit"""

mfdp_keys = [
    "output_file",
    "two_body_interaction",
    "two_body_file_type",
    "Z",
    "N",
    "hbar_omega",
    "Nhw",
    "N_min",
    "N_1max",
    "N_12max",
    "parity",
    "total_2Jz",
    "iham",
    "iclmb",
    "strcm",
    "interaction_type",
    "major",
    "nshll",
    "occupation_string",
    "nsets",
    "min_nesp",
    "nskip",
    "iset1",
    "ki",
    "kf",
    "n_states",
    "gs_energy",
    "iterations_required",
    "igt",
    "irest",
    "nhme",
    "nhw0",
    "nhw_min",
    "nhw_restart",
    "kappa_points",
    "cmin",
    "kappa_restart",
    "kappa_vals",
    "convergence_delta",
    "three_body_interaction",
    "N_123max",
    "eff_charge_p",
    "eff_charge_n",
    "glp",
    "gln",
    "gsp",
    "gsn",
    "saved_pivot",
    "rmemavail"]
"""all variables for mfdp.dat file"""


default_keys = [
    "two_body_file_type",
    "N_min",
    "iham",
    "iclmb",
    "strcm",
    "major",
    "nshll",
    "nsets",
    "min_nesp",
    "nskip",
    "iset1",
    "ki",
    "kf",
    "gs_energy",
    "igt",
    "nhme",
    "nhw_restart",
    "cmin",
    "convergence_delta",
    "eff_charge_p",
    "eff_charge_n",
    "glp",
    "gln",
    "gsp",
    "gsn",
    "mem",
    "tasks_per_node"]
"""all variables that are set by default"""


class Params(object):
    """
    The Params class is pretty much just an argparse.Namespace,
    except you can only have certain variables.

    For example, mfdp.dat files have a field for Z,
    so a Params object for that file type should have an attribute
    self.Z.

    But if you tried to set self.z by accident it should raise an error.

    This prevents accidentally setting the wrong attributes for a file type,
    or too many parameters, or too few.
    """

    def __init__(self, filetype, **kwargs):
        key_map = {
            "MFDP": mfdp_keys,
            "MANUAL INPUT": man_keys,
            "LOCAL_BATCH": local_batch_keys,
            "CEDAR_BATCH": cedar_batch_keys,
            "SUMMIT_BATCH": summit_batch_keys,
            "DEFAULT": default_keys,
            "EMPTY": []}
        self.valid_keys = key_map[filetype]

        # ensure we have the right number of args provided
        if len(kwargs.keys()) < len(self.valid_keys):
            error_message = "Not all necessary parameters were provided to "\
                + filetype + ".\n\nThe missing parameters are:\n"
            for key in self.valid_keys:
                if key not in kwargs.keys():
                    error_message += key+"\n"
            raise ValueError(error_message)
        if len(kwargs.keys()) > len(self.valid_keys):
            error_message = "Too many parameters were provided to "\
                + filetype + ".\n\nThe extra parameters are:\n"
            for key in kwargs.keys():
                if key not in self.valid_keys:
                    error_message += key+"\n"
            raise ValueError(error_message)

        # now ensure all args were valid
        if all(key in self.valid_keys for key in kwargs.keys()):
            # if so, then set self.kwarg = kwarg value
            for key in self.valid_keys:
                setattr(self, key, kwargs[key])
        else:
            error_message = "Not all parameters supplied to " + filetype\
                + " were valid!\n\nThe invalid parameters are:\n"
            for key in kwargs.keys():
                if key not in self.valid_keys:
                    error_message += key+"\n"
            raise ValueError(error_message)

    def param_dict(self):
        """
        Returns a dictionary with the same info contained in the Params object
        """
        pdict = {}
        for key in self.valid_keys:
            pdict[key] = getattr(self, key)
        return pdict


class ManParams(Params):
    """Params subclass for manual input"""
    def __init__(self, **kwargs):
        super(ManParams, self).__init__("MANUAL INPUT", **kwargs)


class LocalBatchParams(Params):
    """Params subclass for local machine batch files"""
    def __init__(self, **kwargs):
        super(LocalBatchParams, self).__init__("LOCAL_BATCH", **kwargs)


class CedarBatchParams(Params):
    """Params subclass for cedar batch files"""
    def __init__(self, **kwargs):
        super(CedarBatchParams, self).__init__("CEDAR_BATCH", **kwargs)


class SummitBatchParams(Params):
    """Params subclass for summit batch files"""
    def __init__(self, **kwargs):
        super(SummitBatchParams, self).__init__("SUMMIT_BATCH", **kwargs)


class MFDPParams(Params):
    """Params subclass for mfdp.dat files"""
    def __init__(self, **kwargs):
        super(MFDPParams, self).__init__("MFDP", **kwargs)


class DefaultParams(Params):
    """Params subclass for default parameters"""
    def __init__(self, **kwargs):
        super(DefaultParams, self).__init__("DEFAULT", **kwargs)


DefaultParamsObj = DefaultParams(
    two_body_file_type=2,
    N_min=0,
    iham=1,
    iclmb=0,
    strcm=10.0,
    major=2,
    nshll=0,
    nsets=3,
    min_nesp=0,
    nskip=0,
    iset1=0,
    ki=1,
    kf=1,
    gs_energy=0.0,
    igt=0,
    nhme=0,
    nhw_restart=-1,
    cmin=0,
    convergence_delta=1.0,
    eff_charge_p=1.0,
    eff_charge_n=0.0,
    glp=1.0,
    gln=0.0,
    gsp=5.586,
    gsn=-3.826,
    mem=0,
    tasks_per_node=48
)
