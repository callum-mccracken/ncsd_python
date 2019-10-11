"""
Holds big long strings to make writing files easier

note: double curly braces tell Python you mean a literal curly brace
      when you're formatting a string.
      "{a}, {b}, {{a}}".format(a=1, b=2) produces the string "1, 2, {a}"
"""

xmgrace_format = """
E\\sx\\N [MeV]
lin lin
{num_spectra_plus_2} {num_states} {num_spectra} {num_plots}              !       num_strings num_states num_specs num_plots
\\S{Z_plus_N}\\N{element}
{interaction_name} 
{axis_labels}
 dotted off rotate
{data}

"""

xmgrace_Nmax_title_format = "N\\smax\\N = {Nmax}"

xmgrace_data_line_format = "{Jx2} {repetition} {parity}   {energy:.4f}\n"

xmgrace_dataset_format = "{title}\n{lines}"

xmgrace_axis_label_line = "{Nmax}h\\h{{-.6}}\\v{{.3}}- \\v{{-.3}}\\xW\\f{{}}\n"
