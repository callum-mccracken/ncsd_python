"""
Contains long strings representing file types, to de-clutter other modules.

PS: if you're reading this online, this looks way better in the actual file

PPS: double curly braces tell Python you mean a literal curly brace
when you're formatting a string.
E.g. "{a}, {b}, {{a}}".format(a=1, b=2) produces the string "1, 2, {a}"
"""

mfdp_format = """  ==== Input file for the shell model code =====
INT filename and output file name:
{two_body_interaction}
{two_body_file_type}     <- TBME file: a b c d  J T  Trel  Hrel  Clmb   V ...
{output_file}
  {Z:2d} {N:2d} {hbar_omega:2.1f}         ! Z, N, hbar*Omega
  {Nhw:2d} {parity:2d} {total_2Jz:2d}           ! Nhw, Parity, Total 2Jz
  {N_min:2d} {N_1max:2d} {N_12max:2d}           ! N_min,N_1max,N_2max
  {iham}  {iclmb}    {strcm}       ! iham(1: H=Trel+V, 2: H=Hrel+V, 3: H=SPE+V), iclmb, strcm
  {interaction_type}                 ! NN=2, NN+3N=3, (NN+3N & NN:IT)=-3
  {major}                 ! major (0 for subshells, 1 for major shells)
  {nshll}                 ! nshll/jorder (norb (major=0) OR orbital order)
Major or subshell occupation restrictions for protons, neutrons and (p+n):
{occupation_string}
Control parameters for multiple-G calculations:
  {nsets}  {min_nesp}  {nskip}  {iset1}         ! nsets (# of G's), min(nesp), nskip, iset1 (see below)
Control parameters for initial and final state vectors:
  {ki}  {kf:2d}  {n_states:2d}          ! ki,kf,nf (Initial, 1st final state, # of final states)
  {gs_energy:.5f}             ! the g.s. energy (0 if not available)
Other control parameters:
  {iterations_required:<5d}               ! Number of lanczos iterations
  {igt}                  ! igt = 1 to calculate sum GT+ and sum GT-, = 0 not to.
Restart option:
  {irest}                  ! irest = 0 for new calc., = 4 to restart.  See code for more options.
  {nhme}                  ! memsave = 1 for a memory saving
Importance truncation:
  {nhw0}  {nhw_min:2d}  {Nhw:2d}  {nhw_restart:2d}      ! nhw0, nhw_min (start IT), nhw_max=nhw, nhw_restart
  {kappa_points}  {cmin}  {kappa_restart}         !kappa_points,cmin=factor*kappa,kappa_restart=-1 for new calc.
  {kappa_vals}  ! kappa(1E-5)
  {convergence_delta}          	     ! convergence delta in keV (not implemented)
{three_body_interaction}
{N_1max:2d} {N_12max:2d} {N_123max:2d}             ! N1B N2B N3B
Effective operators:
    {eff_charge_p:.3f}    {eff_charge_n:.3f}    <- effective charges: ep and en
    {glp:.3f}    {gln:.3f}    {gsp:.3f}   {gsn:.3f}    <- glp,gln,gsp,gsn
{saved_pivot}                    ! pivot_saved
{rmemavail:3f}               ! memory per rank in GB in real(8)
*******************************  End of input  *******************************
"""
"""Format for mfdp.dat file"""

local_batch_format = """#!/bin/bash
cd {run_directory}

potential="{potential}"

iNu="{nucleus_name}"
freq="{hbar_omega}"
suf="{suffix}"
Ngs={Ngs}


{ncsd_path}

for Nmax in {non_IT_Nmax}

do

N=$[$Nmax+$Ngs]

mv mfdp_${{N}}.egv mfdp_${{N}}.egv_${{iNu}}_${{potential}}_Nmax${{Nmax}}.${{freq}}${{suf}}

done
{potential_end_bit}

mv mfd.log mfd.log_{output_file}

python {output_plotter} -f {output_file}
"""
"""Format for local batch file"""

cedar_batch_format = """#!/bin/bash
#SBATCH --account={account}
#SBATCH --nodes={nodes}               # number of 48-cpu nodes
#SBATCH --tasks-per-node={tasks_per_node}      # mpi tasks per node (max 48)
#SBATCH --mem={mem}                 # to use full nodes, set this to zero
#SBATCH --mem-per-cpu={mem_per_core}G     # memory per CPU
#SBATCH --time={time}           # time (DD-HH:MM)
#SBATCH --output={output}

module load python/3.7.0

cd {run_directory}

potential="{potential}"

iNu="{nucleus_name}"
freq="{hbar_omega}"
suf="{suffix}"
Ngs={Ngs}


srun {ncsd_path}

for Nmax in {non_IT_Nmax}

do

N=$[$Nmax+$Ngs]

mv mfdp_${{N}}.egv mfdp_${{N}}.egv_${{iNu}}_${{potential}}_Nmax${{Nmax}}.${{freq}}${{suf}}

done
{potential_end_bit}

mv mfd.log mfd.log_{output_file}

python {output_plotter} -f {output_file}
"""
"""Format for cedar batch file"""

summit_batch_format = """#!/bin/bash

#BSUB -P {account}
#BSUB -W {time}
#BSUB -nnodes {nnodes}
#BSUB -J {output}
#BSUB -eo {output}.%J

cd {run_directory}
date

#module loading needs cuda and for some
#reason also lapack as a separate load
module load cuda
module load netlib-lapack/3.8.0
module load gcc
module load python/3.7.0-anaconda3-5.3.0

#SUMMIT job submission is based on "resource sets"
#so the main thing you submit is "number of resource sets"
#to run efficiently, each resource set has 1 gpu and 7 cpus
#and we run with 6 resource sets per node. This is forced
#by the fact that there are 6 gpus per node and we want each
#to go to a separate resource set.

#full load of node 7 x 6 = 42
export OMP_NUM_THREADS=7

#options:
# -n : number of resource sets TOTAL (nodes x resoure sets per node)
# -a : number of MPI tasks per resource set
# -c : number of CPUs per resource set
# -g : number og GPUs per resource set
# -r : number of resource sets per node
# -l CPU-CPU : optimize CPU memory latency
# -d packed -b rs: pack the resource sets on the node

jsrun -n {resource_sets} -a 1 -c 7 -g 1 -r 6 -l CPU-CPU -d packed -b rs {ncsd_path}

date

potential='{potential}'
iNu='{nucleus_name}'
freq='{hbar_omega}'
suf='{suffix}'
Ngs={Ngs}

for Nmax in {non_IT_Nmax}

do

N=$[$Nmax+$Ngs]

mv mfdp_${{N}}.egv mfdp_${{N}}.egv_${{iNu}}_${{potential}}_Nmax${{Nmax}}.${{freq}}${{suf}}

done

{potential_end_bit}

mv mfd.log mfd.log_{output_file}

python {output_plotter} -f {output_file}
"""
"""Format for summit batch file"""


potential_end_bit_format = """
for Nmax in {IT_Nmax}

do

N=$[$Nmax+$Ngs]

{kappa_rename}
done
"""
"""the ``potential_end_bit`` of batch files"""

kappa_rename_format = """mv mfdp_${{N}}_{kappa_D}.egv mfdp_${{N}}_{kappa_D}.egv_${{iNu}}_${{potential}}_Nmax${{Nmax}}.${{freq}}_IT_kmin{kappa_em}${{suf}}"""
"""format for lines with the kappa renaming scheme"""
