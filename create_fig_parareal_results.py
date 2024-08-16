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
    # Read the namelist
    nml = f90nml.read(config_file_name)
    Tstart = nml['sim_time']['Tstart']
    Tend = nml['sim_time']['Tend']
    do_io_first_coarse = nml['io']['do_io_first_coarse']
    do_io_fine = nml['io']['do_io_fine']
    do_io_cor_coarse = nml['io']['do_io_cor_coarse']
    N_fine = int(nml['num_steps']['N_fine'])
    N_coarse = int(nml['num_steps']['N_coarse'])
    Niter = nml['num_iter']['Niter']
    Nprocs_for_seq = nml['num_proc']['Nprocs_for_seq']
    match_y_axis_range = nml['visualization']['match_y_axis_range']

    #read equential calculation results
    seq_V = np.genfromtxt(seq1_filename)
    
    #pattern of filename
    first_coarse_pattern = "phi_first_coarse_" + "{:02d}".format(Nprocs_for_seq-1) + "_" + "{:02d}".format(Nprocs_for_seq) + "_mpi.dat"
    fine_pattern = "phi_fine_{:02d}_{:02d}_" + "{:02d}".format(Nprocs_for_seq) + "_mpi.dat"
    cor_coarse_pattern = "phi_cor_coarse_{:02d}_{:02d}_" + "{:02d}".format(Nprocs_for_seq) + "_mpi.dat"

    #read parareal calculation results
    rows, coarse_cols, fine_cols = Niter, (N_coarse+1)*Nprocs_for_seq , (N_fine+1)*Nprocs_for_seq
    first_coarse_V = []
    fine_V = [[None for _ in range(fine_cols)] for _ in range(rows)]
    cor_coarse_V = [[None for _ in range(coarse_cols)] for _ in range(rows)]

    if(do_io_first_coarse):
        file_name = first_coarse_pattern
        data = np.genfromtxt(file_name)
        #first_coarse_V.append(data)
        for i in range(len(data)):
            first_coarse_V.append(data[i]) 
    if(do_io_fine):
        for iter in range(1,Niter+1):
            inter_area_counter = 0
            for rank in range(0,Nprocs_for_seq):
                file_name = fine_pattern.format(iter, rank)
                data = np.genfromtxt(file_name)
                for i in range(len(data)):
                    fine_V[iter-1][rank*N_fine + inter_area_counter + i] = data[i]
                inter_area_counter = inter_area_counter + 1
    
                
    if(do_io_cor_coarse):
        for iter in range(1,Niter+1):
            inter_area_counter = 0
            for rank in range(0,Nprocs_for_seq):
                file_name = cor_coarse_pattern.format(iter, rank)
                data = np.genfromtxt(file_name)
                for i in range(len(data)):
                    cor_coarse_V[iter-1][rank*N_coarse + inter_area_counter + i] = data[i]
                inter_area_counter = inter_area_counter + 1

    
    #calculate dt
    seq_dt = (Tend-Tstart) / (N_fine*Nprocs_for_seq)
    coarse_dt = (Tend-Tstart) / (N_coarse*Nprocs_for_seq)
    fine_dt = (Tend-Tstart) / (N_fine*Nprocs_for_seq)

    #Create Time List
    seq_t = []
    first_coarse_t = []
    fine_t = []
    cor_coarse_t = []
    for i in range(len(seq_V)):
        seq_t.append(i*seq_dt)
        
    if(do_io_first_coarse):
        for i in range(len(first_coarse_V)):
            first_coarse_t.append(i*coarse_dt)
    if(do_io_fine):
        for rank in range(0,Nprocs_for_seq):
            for i in range(0,N_fine + 1):
                fine_t.append(rank*N_fine*fine_dt + i*fine_dt)
    if(do_io_cor_coarse):
        for rank in range(0,Nprocs_for_seq):
            for i in range(0,N_coarse + 1):
                cor_coarse_t.append(rank*N_coarse*coarse_dt + i*coarse_dt)
    
    #create figure of potential
    plt.rcParams["font.size"] = 18
    fig1, ax1 = plt.subplots(2,int(Niter/2),figsize = (6*1.6*(Niter/2),12))
    iter = 0
    for row in range(0,2):
        for col in range(0,int(Niter/2)):
            #plot time sequential calculation results
            ax1[row,col].plot(seq_t,seq_V , color='red')
            #plot parareal calculation results
            if(do_io_first_coarse and row == 0 and col == 0):
                ax1[row,col].scatter(first_coarse_t,first_coarse_V , color='orange')
            if(do_io_cor_coarse and (row !=0 or col != 0)):
                ax1[row,col].scatter(cor_coarse_t,cor_coarse_V[iter] , color='orange')
            if(do_io_fine):
                ax1[row,col].scatter(fine_t,fine_V[iter] , color='blue')
            
            if(row == 0 and col == 0 and match_y_axis_range):
                y_axis_limits = ax1[row,col].get_ylim()
            if(match_y_axis_range):
                ax1[row,col].set_ylim(y_axis_limits)
            
            iter = iter + 1
    
    fig1.supxlabel("Time (sec)")
    fig1.supylabel("Potential (V)")
    
    plt.savefig("potential_parareal.png")
    
    #At least memory measures
    first_coarse_V.clear()
    cor_coarse_V.clear()
    first_coarse_t.clear()
    cor_coarse_t.clear()
    
    #Correct data to calculate relative error
    seq_t.pop(0)

    for iter in range(0,Niter):
        for rank in range(Nprocs_for_seq-1, 0, -1):
            fine_V[iter].pop(N_fine*rank + (rank-1))
    
    #calculate relative error
    rel_error = [[0 for _ in range(len(seq_V)-1)] for _ in range(Niter)]
   
    for iter in range (0,Niter):
        for i in range(1,len(seq_V)):
            rel_error[iter][i-1] = abs( (fine_V[iter][i]-seq_V[i]) / seq_V[i] )
    
    #create figure of relative_error
    plt.rcParams["font.size"] = 18
    fig2, ax2 = plt.subplots(2,int(Niter/2),figsize = (6*1.6*(Niter/2),12))
    iter = 0
    for row in range(0,2):
        for col in range(0,int(Niter/2)):
            if(row == 0 and col == 0):
                ax2[row,col].plot(seq_t,rel_error[iter] , color='r')
                ax2[row, col].set_yscale('log')
                y_axis_limits = ax2[row,col].get_ylim()
            else:
                ax2[row,col].plot(seq_t,rel_error[iter] , color='r')
                ax2[row, col].set_yscale('log')
                if(match_y_axis_range):
                    ax2[row,col].set_ylim(y_axis_limits)
                
            iter = iter + 1
    
    fig2.supxlabel("Time (sec)")
    fig2.supylabel("Relative error")

    plt.savefig("relative_error_parareal.png")
    
except Exception as e:
    print(f'Error: {e}')
