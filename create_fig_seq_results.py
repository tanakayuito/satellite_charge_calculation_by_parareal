import f90nml
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import sys
import csv

try:
    #Read inputfile name
    args = sys.argv
    config_file_name = args[1]
    seq1_filename = args[2]
    seq2_filename = args[3]

    # Read the namelist
    nml = f90nml.read(config_file_name)
    do_io_first_coarse = nml['io']['do_io_first_coarse']
    do_io_fine = nml['io']['do_io_fine']
    do_io_cor_coarse = nml['io']['do_io_cor_coarse']
    Tstart = nml['sim_time']['Tstart']
    Tend = nml['sim_time']['Tend']
    N_fine = nml['num_steps']['N_fine']
    N_coarse = nml['num_steps']['N_coarse']
    Niter = nml['num_iter']['Niter']
    Nprocs_for_seq = nml['num_proc']['Nprocs_for_seq']

    #read equential calculation results
    seq1_V = np.genfromtxt(seq1_filename) #results calculated by fine dt
    seq2_V = np.genfromtxt(seq2_filename) #results calculated by coarse dt

    #calculate the length of dt
    seq1_dt = (Tend - Tstart)/ (N_fine*Nprocs_for_seq)
    seq2_dt = (Tend - Tstart)/ (N_coarse*Nprocs_for_seq)
    seq1_dt_list = []
    seq2_dt_list = []
    for i in range(len(seq1_V)):
        seq1_dt_list.append(i*seq1_dt)
    for i in range(len(seq2_V)):
        seq2_dt_list.append(i*seq2_dt)

    #create the figure of potential of seq1
    fig = plt.figure(figsize = (12.8,8))
    plt.rcParams["font.size"] = 18
    plt.plot(seq1_dt_list, seq1_V, color='r')
    plt.xlabel("Time[sec]")
    plt.ylabel("Potential [V]")
    plt.savefig("potential_seq1.png")

    #create the figure of potential of seq2
    fig = plt.figure(figsize = (12.8,8))
    plt.rcParams["font.size"] = 18
    plt.plot(seq2_dt_list, seq2_V, color='r')
    plt.xlabel("Time[sec]")
    plt.ylabel("Potential [V]")
    plt.savefig("potential_seq2.png")

    #correct data for calculating relative error
    seq1_dt_list.pop(0)
    seq2_dt_list.pop(0)

    #calculate reative error
    rel_error = []
    for n in range(1,len(seq2_V)):
        rel_error.append(abs( (seq1_V[n*int(N_fine/N_coarse)]-seq2_V[n]) / seq1_V[n*int(N_fine/N_coarse)]) )

    #create the figure of relative error(seq2 to seq1) 
    fig, ax=plt.subplots(figsize=(12.8,8))
    plt.rcParams["font.size"] = 18
    ax.plot(seq2_dt_list, rel_error, color='r')
    ax.set_xlabel("Time[sec]")
    ax.set_ylabel("Relative error")
    ax.yaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True))
    ax.yaxis.offsetText.set_fontsize(18)
    ax.ticklabel_format(style='sci',axis='y',scilimits=(0,0))

    #if you want to set log scale at yaxis, get out commentout
    ax.set_yscale('log')


    plt.savefig("relative_error_seq2toseq1.png")
    
except Exception as e:
    print(f'Error: {e}')
