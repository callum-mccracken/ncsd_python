"""A module for dealing with reading/writing files for NCSD code."""

from os.path import exists
from . import formats
from . import data_structures
from .data_checker import manual_input_check, check_mfdp_read


class FileManager(object):
    """
    A general class for dealing with files, has a general write function
    that must be implemented for each file type.

    Does not contain a general read function, you have to implement that
    for each different file type if you want one.
    """
    def __init__(self, filetype, filename):
        # the empty case is just good practice, I don't think it's ever used
        # (for much, it might get called in some other default constructors)
        self.params = data_structures.Params("EMPTY")
        self.filename = filename
        if filetype == "MFDP":
            self.valid_keys = data_structures.mfdp_keys
            self.format_string = formats.mfdp_format
        elif filetype == "LOCAL_BATCH":
            self.valid_keys = data_structures.local_batch_keys
            self.format_string = formats.local_batch_format
        elif filetype == "CEDAR_BATCH":
            self.valid_keys = data_structures.cedar_batch_keys
            self.format_string = formats.cedar_batch_format
        elif filetype == "SUMMIT_BATCH":
            self.valid_keys = data_structures.summit_batch_keys
            self.format_string = formats.summit_batch_format
        elif filetype == "DEFAULT":
            self.valid_keys = data_structures.default_keys
            self.format_string = ""
        elif filetype == "EMPTY":
            self.valid_keys = []
            self.format_string = ""

    def param_dict(self):
        """Returns a dictionary containing all file parameters"""
        return self.params.param_dict()

    def write(self):
        """
        Writes ``self.format_string.format(**self.param_dict())``
        to ``self.filename``.
        """
        with open(self.filename, 'w+') as open_file:
            open_file.write(self.format_string.format(**self.param_dict()))


class MFDP(FileManager):
    """
    Subclass of ``FileManager``, for reading / writing ``mfdp.dat`` files.
    """
    def __init__(self, filename="mfdp.dat", params=None):
        super(MFDP, self).__init__("MFDP", filename)
        # if params are given, use those, but if default, read the file
        if params is not None:
            self.params = params
        elif exists(filename):
            self.read()
        elif not exists(filename):
            raise IOError("filename "+filename+" does not exist!")
        else:
            raise IOError("must have either a params object or a filename")

    def read(self):
        """
        Warning:
            This is old, sketchy, and not actually used!
            I just kept it around in case you want to try something similar.
        """
        # open file, grab text
        with open(self.filename, "r") as open_file:
            lines = open_file.readlines()

        # get parameters from each line (brace yourself, this is long...):
        line_num = 0  # for counting lines, but in a weird way
        occupation_string = ""  # we'll need this later too
        for line in lines:
            words = line.split()
            # ignore text lines, except the line that gives saved_pivot!
            if not any(char.isdigit() for char in line):
                if "saved_pivot" in line:
                    saved_pivot = words[0]
            elif line_num == 0:
                two_body_interaction = line
                line_num += 1
            elif line_num == 1:
                two_body_file_type = int(words[0])
                line_num += 1
            elif line_num == 2:
                output_file = line
                line_num += 1
            elif line_num == 3:
                Z = int(words[0])
                N = int(words[1])
                hbar_omega = float(words[2])
                line_num += 1
            elif line_num == 4:
                Nhw = int(words[0])
                parity = int(words[1])
                total_2Jz = int(words[2])
                line_num += 1
            elif line_num == 5:
                N_min = int(words[0])
                N_1max = int(words[1])
                N_12max = int(words[2])
                line_num += 1
            elif line_num == 6:
                iham = int(words[0])
                iclmb = int(words[1])
                strcm = float(words[2])
                line_num += 1
            elif line_num == 7:
                interaction_type = int(words[0])
                line_num += 1
            elif line_num == 8:
                major = int(words[0])
                line_num += 1
            elif line_num == 9:
                nshll = int(words[0])
                line_num += 1
            elif line_num == 10:  # shell diagram thing is of variable length
                if "! N=" in line:
                    occupation_string += line
                else:  # when we're on the next useful line
                    nsets = int(words[0])
                    min_nesp = int(words[1])
                    nskip = int(words[2])
                    iset1 = int(words[3])
                    line_num += 1
            elif line_num == 11:
                ki = int(words[0])
                kf = int(words[1])
                nstates = int(words[2])
                line_num += 1
            elif line_num == 12:
                gs_energy = float(words[0])
                line_num += 1
            elif line_num == 13:
                iterations_required = int(words[0])
                line_num += 1
            elif line_num == 14:
                igt = int(words[0])
                line_num += 1
            elif line_num == 15:
                irest = int(words[0])
                line_num += 1
            elif line_num == 16:
                nhme = int(words[0])
                line_num += 1
            elif line_num == 17:
                nhw0 = int(words[0])
                nhw_min = int(words[1])
                nhw_max = int(words[2])
                nhw_restart = int(words[3])
                line_num += 1
            elif line_num == 18:
                kappa_points = int(words[0])
                cmin = float(words[1])
                kappa_restart = int(words[2])
                line_num += 1
            elif line_num == 19:
                # numbers are words with digits but no letters
                numbers = [
                    word for word in words
                    if (any(char.isdigit() for char in word) and not
                        any(char.isalpha()) for char in word)]
                kappa_vals = " ".join(numbers)
                line_num += 1
            elif line_num == 20:
                convergence_delta = float(words[0])
                line_num += 1
            elif line_num == 21:
                three_body_interaction = line
                line_num += 1
            elif line_num == 22:
                N_1max_verif = int(words[0])
                N_12max_verif = int(words[1])
                N_123max = int(words[2])
                line_num += 1
            elif line_num == 23:
                eff_charge_p = float(words[0])
                eff_charge_n = float(words[1])
                line_num += 1
            elif line_num == 24:
                glp = float(words[0])
                gln = float(words[1])
                gsp = float(words[2])
                gsn = float(words[3])
                line_num += 1
            elif line_num == 25:
                rmemavail = words[0]

        # arrange parameters in a MFDPParams object
        params = data_structures.MFDPParams(
            output_file=output_file,
            two_body_interaction=two_body_interaction,
            two_body_file_type=two_body_file_type,
            Z=Z,
            N=N,
            hbar_omega=hbar_omega,
            Nhw=Nhw,
            N_min=N_min,
            N_1max=N_1max,
            N_12max=N_12max,
            parity=parity,
            total_2Jz=total_2Jz,
            iham=iham,
            iclmb=iclmb,
            strcm=strcm,
            interaction_type=interaction_type,
            major=major,
            nshll=nshll,
            occupation_string=occupation_string,
            nsets=nsets,
            min_nesp=min_nesp,
            nskip=nskip,
            iset1=iset1,
            ki=ki,
            kf=kf,
            n_states=nstates,
            gs_energy=gs_energy,
            iterations_required=iterations_required,
            igt=igt,
            irest=irest,
            nhme=nhme,
            nhw0=nhw0,
            nhw_min=nhw_min,
            nhw_restart=nhw_restart,
            kappa_points=kappa_points,
            cmin=cmin,
            kappa_restart=kappa_restart,
            kappa_vals=kappa_vals,
            convergence_delta=convergence_delta,
            three_body_interaction=three_body_interaction,
            N_123max=N_123max,
            eff_charge_p=eff_charge_p,
            eff_charge_n=eff_charge_n,
            glp=glp,
            gln=gln,
            gsp=gsp,
            gsn=gsn,
            saved_pivot=saved_pivot,
            rmemavail=rmemavail)

        # do a couple quick checks here, then use check_mfdp_read for the rest
        if nhw_max != Nhw:
            raise ValueError("the value of nhw_max must be equal to Nhw")
        if N_1max != N_1max_verif:
            raise ValueError("N_1max from under 3-body file doesn't match")
        if N_12max != N_12max_verif:
            raise ValueError("N_12max from under 3-body file doesn't match")
        check_mfdp_read(params)
        self.params = params


class LocalBatch(FileManager):
    """Subclass of ``FileManager``, for ncsd batch files on local machine."""
    def __init__(self, filename="batch_ncsd", params=None):
        super(LocalBatch, self).__init__("LOCAL_BATCH", filename)
        self.params = params


class CedarBatch(FileManager):
    """Subclass of ``FileManager``, for ncsd batch files on cedar."""
    def __init__(self, filename="batch_ncsd", params=None):
        super(CedarBatch, self).__init__("CEDAR_BATCH", filename)
        self.params = params


class SummitBatch(FileManager):
    """Subclass of ``FileManager``, for ncsd batch files on summit."""
    def __init__(self, filename="batch_ncsd", params=None):
        super(SummitBatch, self).__init__("SUMMIT_BATCH", filename)
        self.params = params


class Defaults(FileManager):
    """Subclass of ``FileManager``,
    for holding default parameters as if they were a file."""
    def __init__(
       self, filename="defaults", params=data_structures.DefaultParamsObj):
        super(Defaults, self).__init__("DEFAULT", filename)
        self.params = params
