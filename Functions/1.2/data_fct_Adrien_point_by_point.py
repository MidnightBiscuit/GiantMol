import numpy as np
from matplotlib import pylab, mlab, pyplot
from pylab import *
from IPython.core.pylabtools import figsize, getfigs
plt = pyplot

import os
from os import listdir
from os.path import isfile, join
import sys

from pathlib import Path
import re
from natsort import natsorted
import timeit

# sys.path.append('/home/adrian/Documents/Programmes/Fortran_Jofre')
# from function_jofre import energy_lost

# sélection des fichiers (Python 3)
import tkinter as tk
from tkinter import filedialog

def load_file_GUI(dir_string):

    root = tk.Tk()
    root.withdraw()

    file_path = filedialog.askopenfilename(initialdir = dir_string,
                                       multiple=True)
    return file_path

def size_data(file_path,row_skip,col_to_read,delim):
    all_size = []
    m = 0
    for k in range(0, len(file_path)):
        temp_size = loadtxt(file_path[k], delimiter=delim,
                            skiprows=row_skip, usecols=col_to_read,
                            unpack=True)[0].size
        all_size.append(temp_size)
        if temp_size > m:
            m = temp_size
    return all_size,m

# import données    /usr/share/myspell/dicts/en_GB.dic
def import_data(file_path,row_skip,col_to_read,delim):

    dico = {}                       # creation dictionnaire
    for i in file_path:              # pour associer à chaque fichier
        dico[i] = []                # un nombre de colonnes dépendant de col_to_read
# dico de la forme
# dico = {'my_file1.csv': array([1,2,3]),
#          'my_file2.csv': array([2,4,6]),
#          'my_file3.csv': array([5,10,15])}
# récup valeur dans dico : dico.get('my_file1.csv')

    for k in range(0, len(file_path)):
        dico[file_path[k]] = loadtxt(file_path[k], delimiter=delim,encoding ='utf-8',
                                      skiprows=row_skip, usecols=col_to_read,
                                      unpack=True)
    return dico

def convert_dico_to_var(dico):

    file_path = list(dico.keys()) # all the keys from dico

    # création variables courantes
    n = len(dico)  # nombre de fichiers
    m = len(dico.get(file_path[0])[0])  # longueur colones
    p = len(dico.get(file_path[0])[1:])  # nombre de col : temps + canaux oscillo
    shape = (n, m)
    TP = zeros(shape)  # time
    CH = zeros((n, p, m))  # les canaux de l'oscillo (p-1)  CH[file:channel:time]

    for k in range(0, n):
        TP[k, :] = dico.get(file_path[k])[0, :]
        for l in range(1, p+1):
            CH[k, l-1, :] = dico.get(file_path[k])[l, :]

    return TP,CH

##### FUNCTION_JOFRE.IPY #####

def load_T_and_PM_simu(str_load):
    data = loadtxt('{}.dat'.format(str_load),comments='%');
    return data[:,0],data[:,1:4], data[:,4:7], data[:,7] # tt, T_CM, T_aux, PM

def load_xyz_init_bin_DP(str_load):
    aux_info = loadtxt(str_load+'.info',comments='%');
    n_ions1 = int(aux_info[0])
    n_ions2 = int(aux_info[1])
    n_ions = n_ions1 + n_ions2
#    print(n_ions1, n_ions2)
    imax = 12 * n_ions

#    print('extra_time', aux_info[19])
    
    fid = open(str_load+'.bin', 'rb')
    junk = fromfile(fid, int32,1)        # Read record start tag
    aux  = fromfile(fid, int32, 3);
    junk = fromfile(fid, int32,1)        # Read record stop tag

    junk = fromfile(fid, int32,1)        # Read record start tag
    aux  = fromfile(fid, float64, imax);
    junk = fromfile(fid, int32,1)        # Read record stop tag

    fid.close
    aux = reshape(aux, (12, n_ions),order='F')
    r_LC = 1.e3*aux[0:3,0:n_ions1]
    v_LC = aux[3:6,0:n_ions1]

    #xyz = 1.e3*aux[0:3,:]
    a_LC = aux[6:9,0:n_ions1]
    v_LC_avg = aux[9:12,0:n_ions1]

    return r_LC,v_LC,a_LC,v_LC_avg

def plot_XYZ(file_name,fig_name='2',fig_title='XYZ'):
    r_LC,v_LC,a_LC,v_LC_avg = load_xyz_init_bin_DP(file_name)
    figure(fig_name); clf()
    title(fig_title)
    subplot(211,aspect=1.0)
    plot(r_LC [0,:],r_LC [1,:],'8',color='xkcd:purplish blue')
    xlabel('x[mm]')
    ylabel('y[mm]')
    grid()

    # subplot(212,aspect=1.0)
    subplot(212)
    plot(r_LC [2,:],r_LC [0,:],'8',color='xkcd:purplish blue')
    xlabel('z[mm]')
    ylabel('x[mm]')
    grid()
    
# def plot_T_and_PM_Init_Inje_Evol(file_dir,file1, file2, file3,fig_name='3'):

def plot_T_and_PM_Init_Inje_Evol(file_dir2,file_name,flag_plot,fig_name,**kwargs):
    
    # ~ xlim1 = (-0.1,6)
    # ~ ylim1 = (0.5*1e-3,5e3)
    # ~ ylim2 = (-2,120)
            
    xlim1 = kwargs.get('xlim1', (-0.1,6))
    ylim1 = kwargs.get('ylim1', (0.5*1e-3,2e4))
    ylim2 = kwargs.get('ylim2', (-2,50))
    
    i_aux = file_name.find('_N')
    # file0 = 'Temp_3D_Harmo'+ file_name[i_aux:]
    file1 = 'SimuType0'    + file_name[i_aux:] ########## file1 = 'SimuType0'    + file_name[i_aux:17+36]
    file2 = 'SimuType4_01' + file_name[i_aux:]
    file3 = 'SimuType2_01' + file_name[i_aux:]

    # tt0, T_CM0, T_aux0, PM0 = load_T_and_PM_simu(file_dir2+file0)
    tt1, T_CM1, T_aux1, PM1 = load_T_and_PM_simu(file_dir2+'Temp_'+file1)
    tt2, T_CM2, T_aux2, PM2 = load_T_and_PM_simu(file_dir2+'Temp_'+file2+'50eV')
    tt3, T_CM3, T_aux3, PM3 = load_T_and_PM_simu(file_dir2+'Temp_'+file3+'50eV')

    aux = mean(PM1[-100:])
    PM_variation = ( aux - mean(PM3[-100:]) ) / aux
    T_variation  = mean(T_aux3[-100:,0]) + mean(T_aux3[-100:,1]) + mean(T_aux3[-100:,2])
            

    # Auxiliary arrays:
    t_aux1 = array([tt2[ 0],tt2[ 0]])
    t_aux2 = array([tt2[-1],tt2[-1]])
    y1_aux = array([1.0e-3 ,1.0e-1 ])
#     y2_aux = array([0 ,20 ])
    y2_aux = array([0 ,50 ])

    tt    = concatenate( (   tt1,   tt2,   tt3) )
    T_CM  = concatenate( (T_CM1,T_CM2,T_CM3) )
    T_aux = concatenate( (T_aux1,T_aux2,T_aux3) )
    PM    = concatenate( (PM1,PM2,PM3) )
    SNR = np.abs( aux - mean(PM3[-100:]) )/np.sqrt(aux) # Sig-to-Noi ratio considering Poisson noise
    
    if flag_plot == 1 :
        #fig_name = file_name[-9:]
        fig = figure(fig_name); clf()
        ax1 = subplot(211)
        semilogy(tt*1.e3,T_aux[:,0], label='Tx')
        semilogy(tt*1.e3,T_aux[:,1], label='Ty')
        semilogy(tt*1.e3,T_aux[:,2], label='Tz')
        semilogy(t_aux1*1.e3,y1_aux,'r')
        semilogy(t_aux2*1.e3,y1_aux,'r')
        ax1.grid()

        legend()
        # ~ xlabel('time[ms]')
        # ~ ylabel('T[K]')
        plt.setp(ax1.get_xticklabels(),visible=False)

        ax2 = subplot(212,sharex=ax1)
        plot(tt*1.e3,PM[:])
        plot(t_aux1*1.e3,y2_aux,'r')
        plot(t_aux2*1.e3,y2_aux,'r')
        ax2.grid()
        
        xlabel('time[ms]')
        ylabel('Counts')
                               
        ax1.set_xlim(xlim1)
        ax1.set_ylim(ylim1)
        ax2.set_ylim(ylim2)
        plt.tight_layout()
        subplots_adjust(hspace=0.015)
        
    temps = [tt,t_aux1, t_aux2]
    fluo = [PM, PM_variation, SNR]
    temperature = [T_aux, T_variation, T_CM]
        
    return temps, temperature, fluo
    
def plot_T_and_PM_InitQ_Inje_Evol(file_dir2,file_name,flag_plot,fig_name,**kwargs):
    
    # ~ xlim1 = (-0.1,6)
    # ~ ylim1 = (0.5*1e-3,5e3)
    # ~ ylim2 = (-2,120)
            
    xlim1 = kwargs.get('xlim1', (-0.1,6))
    ylim1 = kwargs.get('ylim1', (0.5*1e-3,2e4))
    ylim2 = kwargs.get('ylim2', (-2,50))
    
    i_aux = file_name.find('_N')
    # file0 = 'Temp_3D_Harmo'+ file_name[i_aux:]
    file1 = 'SimuTypeQ'    + file_name[i_aux:] ########## file1 = 'SimuType0'    + file_name[i_aux:17+36]
    file2 = 'SimuType4_01' + file_name[i_aux:]
    file3 = 'SimuType2_01' + file_name[i_aux:]

    # tt0, T_CM0, T_aux0, PM0 = load_T_and_PM_simu(file_dir2+file0)
    tt1, T_CM1, T_aux1, PM1 = load_T_and_PM_simu(file_dir2+'Temp_'+file1)
    tt2, T_CM2, T_aux2, PM2 = load_T_and_PM_simu(file_dir2+'Temp_'+file2+'50eV')
    tt3, T_CM3, T_aux3, PM3 = load_T_and_PM_simu(file_dir2+'Temp_'+file3+'50eV')

    aux = mean(PM1[-100:])
    PM_variation = ( aux - mean(PM3[-100:]) ) / aux
    T_variation  = mean(T_aux3[-100:,0]) + mean(T_aux3[-100:,1]) + mean(T_aux3[-100:,2])
            

    # Auxiliary arrays:
    t_aux1 = array([tt2[ 0],tt2[ 0]])
    t_aux2 = array([tt2[-1],tt2[-1]])
    y1_aux = array([1.0e-3 ,1.0e-1 ])
#     y2_aux = array([0 ,20 ])
    y2_aux = array([0 ,50 ])

    tt    = concatenate( (   tt1,   tt2,   tt3) )
    T_CM  = concatenate( (T_CM1,T_CM2,T_CM3) )
    T_aux = concatenate( (T_aux1,T_aux2,T_aux3) )
    PM    = concatenate( (PM1,PM2,PM3) )
    SNR = np.abs( aux - mean(PM3[-100:]) )/np.sqrt(aux) # Sig-to-Noi ratio considering Poisson noise
    
    if flag_plot == 1 :
        #fig_name = file_name[-9:]
        fig = figure(fig_name); clf()
        ax1 = subplot(211)
        semilogy(tt*1.e3,T_aux[:,0], label='Tx')
        semilogy(tt*1.e3,T_aux[:,1], label='Ty')
        semilogy(tt*1.e3,T_aux[:,2], label='Tz')
        semilogy(t_aux1*1.e3,y1_aux,'r')
        semilogy(t_aux2*1.e3,y1_aux,'r')
        ax1.grid()

        legend()
        # ~ xlabel('time[ms]')
        # ~ ylabel('T[K]')
        plt.setp(ax1.get_xticklabels(),visible=False)

        ax2 = subplot(212,sharex=ax1)
        plot(tt*1.e3,PM[:])
        plot(t_aux1*1.e3,y2_aux,'r')
        plot(t_aux2*1.e3,y2_aux,'r')
        ax2.grid()
        
        xlabel('time[ms]')
        ylabel('Counts')
                               
        ax1.set_xlim(xlim1)
        ax1.set_ylim(ylim1)
        ax2.set_ylim(ylim2)
        plt.tight_layout()
        subplots_adjust(hspace=0.015)
        
    temps = [tt,t_aux1, t_aux2]
    fluo = [PM, PM_variation, SNR]
    temperature = [T_aux, T_variation, T_CM]
        
    return temps, temperature, fluo

def plot_T_and_PM_Init(file_dir2,file_name,**kwargs):
    
    # ~ xlim1 = (-0.1,6)
    # ~ ylim1 = (0.5*1e-3,5e3)
    # ~ ylim2 = (-2,120)
            
    xlim1 = kwargs.get('xlim1', (-0.1,6))
    ylim1 = kwargs.get('ylim1', (0.5*1e-3,2e4))
    ylim2 = kwargs.get('ylim2', (-2,50))
    
    i_aux = file_name.find('_N')
    # file0 = 'Temp_3D_Harmo'+ file_name[i_aux:]
    file1 = 'SimuType0'    + file_name[i_aux:] ########## file1 = 'SimuType0'    + file_name[i_aux:17+36]

    # tt0, T_CM0, T_aux0, PM0 = load_T_and_PM_simu(file_dir2+file0)
    tt1, T_CM1, T_aux1, PM1 = load_T_and_PM_simu(file_dir2+'Temp_'+file1)
            

    # Auxiliary arrays:
    t_aux1 = array([tt1[-1],tt1[-1]])
    y1_aux = array([1.0e-3 ,1.0e-1 ])
        
    return tt1, T_CM1, T_aux1, PM1

def plot_T_and_PM_Quen(file_dir2,file_name,**kwargs):
    
    # ~ xlim1 = (-0.1,6)
    # ~ ylim1 = (0.5*1e-3,5e3)
    # ~ ylim2 = (-2,120)
            
    xlim1 = kwargs.get('xlim1', (-0.1,6))
    ylim1 = kwargs.get('ylim1', (0.5*1e-3,2e4))
    ylim2 = kwargs.get('ylim2', (-2,50))
    
    i_aux = file_name.find('_N')
    # file0 = 'Temp_3D_Harmo'+ file_name[i_aux:]
    file1 = 'SimuTypeQ'    + file_name[i_aux:] ########## file1 = 'SimuType0'    + file_name[i_aux:17+36]

    # tt0, T_CM0, T_aux0, PM0 = load_T_and_PM_simu(file_dir2+file0)
    tt1, T_CM1, T_aux1, PM1 = load_T_and_PM_simu(file_dir2+'Temp_'+file1)
            

    # Auxiliary arrays:
    t_aux1 = array([tt1[-1],tt1[-1]])
    y1_aux = array([1.0e-3 ,1.0e-1 ])
        
    return tt1, T_CM1, T_aux1, PM1

def plot_T_and_PM_Inje(file_dir2,file_name,**kwargs):
   
    i_aux = file_name.find('_N')
    file2 = 'SimuType4_01' + file_name[i_aux:]

    tt2, T_CM2, T_aux2, PM2 = load_T_and_PM_simu(file_dir2+'/Temp_'+file2+'50eV')

    T_variation  = mean(T_aux2[-25:,0]) + mean(T_aux2[-25:,1]) + mean(T_aux2[-25:,2])

    # Auxiliary arrays:
    t_aux1 = array([tt2[ 0],tt2[ 0]])
    t_aux2 = array([tt2[-1],tt2[-1]])
    y1_aux = array([1.0e-3 ,1.0e-1 ])
#     y2_aux = array([0 ,20 ])
    y2_aux = array([0 ,50 ])
    
    return tt2, T_CM2, T_aux2, PM2

def plot_T_and_PM_Evol(file_dir2,file_name,**kwargs):

    i_aux = file_name.find('_N')
    file3 = 'SimuType2_01' + file_name[i_aux:]

    tt3, T_CM3, T_aux3, PM3 = load_T_and_PM_simu(file_dir2+'Temp_'+file3+'50eV')

    T_variation  = mean(T_aux3[-100:,0]) + mean(T_aux3[-100:,1]) + mean(T_aux3[-100:,2])
        
    return tt3, T_CM3, T_aux3, PM3

def plot_T_and_PM_Init_RFrelax_AfterCooling_Evol(file_dir2,file_name,flag_plot,fig_name,**kwargs):
    
    # ~ xlim1 = (-0.1,6)
    # ~ ylim1 = (0.5*1e-3,5e3)
    # ~ ylim2 = (-2,120)
            
    xlim1 = kwargs.get('xlim1', (-0.1,6))
    ylim1 = kwargs.get('ylim1', (0.5*1e-3,2e4))
    ylim2 = kwargs.get('ylim2', (-2,50))
    
    i_aux = file_name.find('_N')
    file1 = 'SimuType0'    + file_name[i_aux:]
    file2 = 'SimuType4_01' + file_name[i_aux:]
    file3 = 'SimuType2_01' + file_name[i_aux:]
    file4 = 'SimuType6_01' + file_name[i_aux:]

    tt1, T_CM1, T_aux1, PM1 = load_T_and_PM_simu(file_dir2+'Temp_'+file1)
    tt2, T_CM2, T_aux2, PM2 = load_T_and_PM_simu(file_dir2+'Temp_'+file2+'50eV')
    tt3, T_CM3, T_aux3, PM3 = load_T_and_PM_simu(file_dir2+'Temp_'+file3+'50eV')
    try:
        tt4, T_CM4, T_aux4, PM4 = load_T_and_PM_simu(file_dir2+'Temp_'+file4+'50eV')        
    except:
        tt4, T_CM4, T_aux4, PM4 = load_T_and_PM_simu(file_dir2+'Temp_'+file4)

    T_variation  = mean(T_aux2[-100:,0])+mean(T_aux2[-100:,1])+mean(T_aux2[-100:,2]) - mean(T_aux1[-100:,0])+mean(T_aux1[-100:,1])+mean(T_aux1[-100:,2])
    aux = mean(PM1[-100:])
    PM_variation = ( aux - mean(PM3[-100:]) ) / aux
    SNR = np.abs( aux - mean(PM3[-100:]) )/np.sqrt(aux) # Sig-to-Noi ratio considering Poisson noise

    # Auxiliary arrays:
    t_aux1 = array([tt4[ 0],tt4[ 0]])
    t_aux2 = array([tt4[-1],tt4[-1]])
    y1_aux = array([1.0e-4 ,3.0e4 ])
#     y2_aux = array([0 ,20 ])
    y2_aux = array([0 ,50 ])

    tt    = concatenate( (   tt1,   tt4) )
    T_CM  = concatenate( (T_CM1, T_CM4) )
    T_aux = concatenate( (T_aux1,T_aux4) )
    PM    = concatenate( (PM1,PM4) )
    
    if flag_plot == 1 :
        #fig_name = file_name[-9:]
        fig = figure(fig_name); clf()
        ax1 = subplot(111)
        semilogy(tt*1.e3,T_aux[:,0], label='Tx')
        semilogy(tt*1.e3,T_aux[:,1], label='Ty')
        semilogy(tt*1.e3,T_aux[:,2], label='Tz')
        semilogy(t_aux1*1.e3,y1_aux,'r')
        semilogy(t_aux2*1.e3,y1_aux,'r')
        ax1.grid()
        # ax1.set_xticks([1e-4,1e-2,1,1e2,1e4])

        legend()
        xlabel('time[ms]')
        ylabel('T[K]')
        
        fig.set_size_inches(11.69, 8.27)
        plt.tight_layout()
        
    temps = [tt,t_aux1, t_aux2]
    fluo = [PM, PM_variation, SNR]
    temperature = [T_aux, T_variation, T_CM]
        
    return temps, temperature, fluo

def plot_T_and_PM_Init_RFrelax_AfterInj_Evol(file_dir2,file_name,flag_plot,fig_name,**kwargs):
    
    # ~ xlim1 = (-0.1,6)
    # ~ ylim1 = (0.5*1e-3,5e3)
    # ~ ylim2 = (-2,120)
            
    xlim1 = kwargs.get('xlim1', (-0.1,6))
    ylim1 = kwargs.get('ylim1', (0.2*1e-3,2e4))
    ylim3 = kwargs.get('ylim2', (-2,85))
    
    i_aux = file_name.find('_N')
    file1 = 'SimuType0'    + file_name[i_aux:] ########## file1 = 'SimuType0'    + file_name[i_aux:17+36]
    file2 = 'SimuType4_01' + file_name[i_aux:]
    file3 = 'SimuType2_01' + file_name[i_aux:]
    file4 = 'SimuType6_01' + file_name[i_aux:]

    tt1, T_CM1, T_aux1, PM1 = load_T_and_PM_simu(file_dir2+'Temp_'+file1)
    tt2, T_CM2, T_aux2, PM2 = load_T_and_PM_simu(file_dir2+'Temp_'+file2+'50eV')
    tt3, T_CM3, T_aux3, PM3 = load_T_and_PM_simu(file_dir2+'Temp_'+file3+'50eV')
    try:
        tt4, T_CM4, T_aux4, PM4 = load_T_and_PM_simu(file_dir2+'Temp_'+file4+'50eV')        
    except:
        tt4, T_CM4, T_aux4, PM4 = load_T_and_PM_simu(file_dir2+'Temp_'+file4)
    
    aux = mean(PM1[-100:])
    PM_variation = ( aux - mean(PM3[-100:]) ) / aux
    T_variation  = mean(T_aux3[-100:,0]) + mean(T_aux3[-100:,1]) + mean(T_aux3[-100:,2])
    SNR = np.abs( aux - mean(PM3[-100:]) )/np.sqrt(aux) # Sig-to-Noi ratio considering Poisson noise

    # Auxiliary arrays:
    t_aux1 = array([tt2[ 0],tt2[ 0]])
    t_aux2 = array([tt2[-1],tt2[-1]])
    y1_aux = array([1.0e-3 ,1.0 ])
#     y2_aux = array([0 ,20 ])
    y2_aux = array([0 ,50 ])

    tt    = concatenate( (   tt1,   tt2,   tt3) )
    T_aux = concatenate( (T_aux1,T_aux2,T_aux3) )
    T_CM  = concatenate( (T_CM1,T_CM2,T_CM3) )
    PM    = concatenate( (PM1,PM2,PM3) )
    
    tt_relax    = concatenate( (   tt1,   tt2,   tt4) )
    T_aux_relax = concatenate( (T_aux1,T_aux2,T_aux4) )
    PM_relax    = concatenate( (PM1,PM2,PM4) )
    
    if flag_plot == 1 :
        #fig_name = file_name[-9:]
        fig = figure(fig_name); clf()
        ax1 = subplot(311)
        semilogy(tt*1.e3,T_aux[:,0], label='Tx')
        semilogy(tt*1.e3,T_aux[:,1], label='Ty')
        semilogy(tt*1.e3,T_aux[:,2], label='Tz')
        semilogy(t_aux1*1.e3,y1_aux,'r')
        semilogy(t_aux2*1.e3,y1_aux,'r')
        ax1.set_yticks([1e-4,1e-2,1,1e2,1e4])
        ax1.grid()
        # annotate('Laser ON', xy=(0.5,350), xycoords='data',
            # size=24, ha='left', va='top', color='xkcd:azul',
            # bbox=dict(boxstyle='round', fc='white',edgecolor='xkcd:azul'))
        

        legend(title='Laser ON',fontsize=18)
        # ~ xlabel('time[ms]')
        # ~ ylabel('T[K]')
        plt.setp(ax1.get_xticklabels(),visible=False)

        ax2 = subplot(312,sharex=ax1,sharey=ax1)
        semilogy(tt_relax*1.e3,T_aux_relax[:,0], label='Tx',nonposy='mask')
        semilogy(tt_relax*1.e3,T_aux_relax[:,1], label='Ty',nonposy='mask')
        semilogy(tt_relax*1.e3,T_aux_relax[:,2], label='Tz',nonposy='mask')
        semilogy(t_aux1*1.e3,y1_aux,'r')
        semilogy(t_aux2*1.e3,y1_aux,'r')
        ax2.set_yticks([1e-4,1e-2,1,1e2,1e4])
        ax2.grid()
        # annotate('laser off après injection', xy=(0.5,350), xycoords='data',
            # size=24, ha='left', va='top', color='xkcd:azul',
            # bbox=dict(boxstyle='round', fc='white',edgecolor='xkcd:azul'))
        legend(title='Laser OFF après injection',fontsize=18)    
        plt.setp(ax2.get_xticklabels(),visible=False)

        ax3 = subplot(313,sharex=ax1)
        plot(tt*1.e3,PM[:])
        plot(t_aux1*1.e3,y2_aux,'r')
        plot(t_aux2*1.e3,y2_aux,'r')
        ax3.grid()
        
        xlabel('time[ms]')
        ylabel('Counts')
                               
        ax1.set_xlim(xlim1)
        ax1.set_ylim(ylim1)
        ax3.set_ylim(ylim3)
        plt.tight_layout()
        
        fig.set_size_inches(11.69, 8.27)
        subplots_adjust(hspace=0.015)
        
        temps = [tt,t_aux1, t_aux2]
        fluo = [PM, PM_variation, SNR]
        temperature = [T_aux, T_variation, T_CM]
        
    return temps, temperature, fluo

def plot_T_and_PM_InitQ_Evol_AfterCool(data, flag_plot, fig_name, **kwargs):
    onlyfiles = sort([f for f in listdir(data) if isfile(join(data, f)) and "SimuType" in f and ".dat" in f])
    onlyfiles_Lan = sort([f for f in listdir(data) if isfile(join(data, f)) and "Harmo" in f])

    xlim1 = kwargs.get('xlim1', (-0.1, 6))
    ylim1 = kwargs.get('ylim1', (0.5 * 1e-3, 2e4))
    ylim2 = kwargs.get('ylim2', (-2, 50))

    N_ions, j_save, dt_j_save_next, eta, Temp, save_T = \
        load_Temp_init_bin_Lan(data + os.path.sep + onlyfiles_Lan[0].strip('.bin'), 0)
    tt1, T_CM1, T_aux1, PM1 = load_T_and_PM_simu(data + os.path.sep + onlyfiles[1].strip('.dat'))
    tt3, T_CM3, T_aux3, PM3 = load_T_and_PM_simu(data + os.path.sep + onlyfiles[0].strip('.dat'))

    tt0 = save_T[:, 0]
    T_CM0 = save_T[:, 1:4]
    T_aux0 = save_T[:, 4:]
    tt = concatenate((tt0, tt1, tt3))
    T_CM = concatenate((T_CM0, T_CM1, T_CM3))
    T_aux = concatenate((T_aux0, T_aux1, T_aux3))

    if flag_plot == 1:
        # fig_name = file_name[-9:]
        fig = figure(fig_name);
        clf()
        ax1 = subplot(211)
        semilogy(tt * 1.e3, T_aux[:, 0], label='Tx')
        semilogy(tt * 1.e3, T_aux[:, 1], label='Ty')
        semilogy(tt * 1.e3, T_aux[:, 2], label='Tz')
        ax1.grid()

        legend()
        # ~ xlabel('time[ms]')
        # ~ ylabel('T[K]')
        plt.setp(ax1.get_xticklabels(), visible=False)

        ax2 = subplot(212, sharex=ax1)
        plot(tt * 1.e3, T_aux[:, 0], label='Tx')
        plot(tt * 1.e3, T_aux[:, 1], label='Tx')
        plot(tt * 1.e3, T_aux[:, 2], label='Tx')
        ax2.grid()

        xlabel('time[ms]')
        ylabel('T [k]')

        ax1.set_xlim(xlim1)
        ax1.set_ylim(ylim1)
        ax2.set_ylim(ylim2)
        plt.tight_layout()
        subplots_adjust(hspace=0.015)

    return tt, T_aux, T_CM


def find_PM_variation_FinalT(file_dir2,file_name):
    i_aux = file_name.find('_N')+1
    file1 = 'SimuTypeQ_{}'.format(file_name[i_aux:].strip('50eV.dat'))    # for N=1024
    file3 = 'SimuType2_01_{}50eV'.format(file_name[i_aux:].strip('50eV.dat'))
    
    tt1, T_CM1, T_aux1, PM1 = load_T_and_PM_simu(file_dir2+'/Temp_'+file1)
    tt3, T_CM3, T_aux3, PM3 = load_T_and_PM_simu(file_dir2+'/Temp_'+file3)
    
    aux = mean(PM1[-100:])
    PM_variation = ( aux - mean(PM3[-100:]) ) / aux
    T_variation  = mean(T_aux3[-100:,0]) + mean(T_aux3[-100:,1]) + mean(T_aux3[-100:,2])
    
    SNR = np.abs( aux - mean(PM3[-100:]) )/np.sqrt(aux) # Sig-to-Noi ratio considering Poisson noise
    
    return PM_variation, T_variation, SNR

# Look at the lost of energy:
def energy_lost(file_dir2,file_name):
    
    data = loadtxt('{}/{}.dat'.format(file_dir2,file_name))

    r_0 = data[0:3]; v_0 = data[3: 6]; # Position and velocity at the start
    r_1 = data[6:9]; v_1 = data[9:12]; # Position and velocity at the end
    t_c = data[12] ; # Time of crossing

    Ec0 = sum(v_0*v_0) ; Ec1 = sum(v_1*v_1)
    delta_Ec = (Ec1 - Ec0)
    return r_0, r_1, delta_Ec, delta_Ec / Ec1, t_c
    
def load_trj(str_load):
   
    fid = open(str_load+'.bin', 'rb')

    junk = fromfile(fid, int32,1)        # Read record start tag
    aux  = fromfile(fid, int32,1)
    junk = fromfile(fid, int32,1)        # Read record stop tag
    jmax = aux[0]

    junk = fromfile(fid, int32,1)        # Read record start tag
    aux  = fromfile(fid, int32,1)
    junk = fromfile(fid, int32,1)        # Read record stop tag
    N_ions1 = aux[0]

    junk = fromfile(fid, int32,1)        # Read record start tag
    aux  = fromfile(fid, int32,1)
    junk = fromfile(fid, int32,1)        # Read record stop tag
    N_ions2 = aux[0]
 
    junk = fromfile(fid, int32  ,1   )        # Read record start tag
    t    = fromfile(fid, float64,jmax)
    junk = fromfile(fid, int32  ,1   )        # Read record stop tag

    junk = fromfile(fid, int32  ,1   )        # Read record start tag
    r_x  = fromfile(fid, float64,jmax)
    junk = fromfile(fid, int32  ,1   )        # Read record stop tag

    junk = fromfile(fid, int32  ,1   )        # Read record start tag
    r_y  = fromfile(fid, float64,jmax)
    junk = fromfile(fid, int32  ,1   )        # Read record stop tag

    junk = fromfile(fid, int32  ,1   )        # Read record start tag
    r_z  = fromfile(fid, float64,jmax)
    junk = fromfile(fid, int32  ,1   )        # Read record stop tag

    junk = fromfile(fid, int32  ,1   )        # Read record start tag
    v_x  = fromfile(fid, float64,jmax)
    junk = fromfile(fid, int32  ,1   )        # Read record stop tag

    junk = fromfile(fid, int32  ,1   )        # Read record start tag
    v_y  = fromfile(fid, float64,jmax)
    junk = fromfile(fid, int32  ,1   )        # Read record stop tag

    junk = fromfile(fid, int32  ,1   )        # Read record start tag
    v_z  = fromfile(fid, float64,jmax)
    junk = fromfile(fid, int32  ,1   )        # Read record stop tag

    junk = fromfile(fid, int32  ,1   )        # Read record start tag
    a_x  = fromfile(fid, float64,jmax)
    junk = fromfile(fid, int32  ,1   )        # Read record stop tag

    junk = fromfile(fid, int32  ,1   )        # Read record start tag
    a_y  = fromfile(fid, float64,jmax)
    junk = fromfile(fid, int32  ,1   )        # Read record stop tag

    junk = fromfile(fid, int32  ,1   )        # Read record start tag
    a_z  = fromfile(fid, float64,jmax)
    junk = fromfile(fid, int32  ,1   )        # Read record stop tag

    fid.close

    return t,r_x, r_y, r_z, v_x, v_y, v_z, a_x, a_y, a_z


def load_cloud_trj(str_load):
    fid = open(str_load + '.bin', 'rb')

    junk0 = fromfile(fid, int32, 1)  # Read record start tag
    aux = fromfile(fid, int32, 1)
    junk1 = fromfile(fid, int32, 1)  # Read record stop tag
    jmax = aux[0]

    junk2 = fromfile(fid, int32, 1)  # Read record start tag
    aux = fromfile(fid, int32, 1)
    junk3 = fromfile(fid, int32, 1)  # Read record stop tag
    N_ions = aux[0]

    junk4 = fromfile(fid, int32, 1)  # Read record start tag
    r = fromfile(fid, single, int(junk4[0] / 4))  # jmax*N_ions*3
    junk5 = fromfile(fid, int32, 1)  # Read record stop tag

    fid.close

    return jmax, N_ions, r
    
##### FORTRAN ANALYSIS #####
###   Fonctions Adrien   ###

def load_gui(filter_nocomplete):

    # sélection des fichiers (Python 3)
    # SELECTIONNER UN FICHIER DE DONNÉES AU FOND DE L'ARBORESCENCE
    # Temp_SimuType0_N01024_Vrf0064_Udc0.5000D+00V_D1.0_S1.0RFG.dat

#     root = tk.Tk()
#     root.withdraw()

#     # Sélectionner un fichier du dossier pour récupérer le répertoire et les noms
#     file_path = filedialog.askopenfilename(initialdir = '/home/adrian/Documents/Programmes/Fortran_Jofre',
#                                            multiple=None)
    
    file_path = load_file_GUI('/home/adrian/Documents/Programmes/Fortran_Jofre')[0]
    
    dir_path = Path(file_path).parent.absolute()
    work_rep = file_path[:len(str(dir_path))]
    filename = file_path[len(str(dir_path))+1:]

    print('{}{}'.format('> Répertoire : ',work_rep))
    print('{}{}'.format('> Filename : ',filename))
    
    ## Recover all subdirectories hosting simulation data

    # find all '/' in path
    myslashpos = [m.start() for m in re.finditer('/', work_rep)]
    # print(myslashpos[6]) # choose the right level

    if 'home' in file_path:
        if 'Hobitton' or 'Rivendel' in file_path :
            slashcond = 5 # the slash number just before runs (after date)
        else:
            slashcond = 6 # the slash number just before runs (after date)
    else:
        slashcond = -2

    print('> myslashpos |',myslashpos)
    print('> slashcond |',slashcond)

    # Getting a naturally sorted list of all subdirectories in some directory
    # throw all the intermediary subdirectories that are not hosting data
    all_subdir = [x[0] for x in os.walk(work_rep[0:myslashpos[slashcond]])]
    all_subdir = natsorted(all_subdir)

    # keep only the long path (by remove short ones)
    all_subdir = [i for i in all_subdir if len(i) >= len(str(dir_path))-1]

    if filter_nocomplete == 1:
        # Name of all points (conditions) 
        all_points = [point[myslashpos[slashcond] + 1:myslashpos[slashcond + 1]] for point in all_subdir]
        all_points = list(dict.fromkeys(all_points))
        # how much repetition for first point : usually the first point is complete
        # i'm not going to plot the data while not even the first point is complete ...
        num_rep = len([k for k in all_subdir if all_points[0] in k])
        
        # determine which points are not complete
        # the way the code is done usually leads to
        # only one uncomplete point but in any case
        # I handle the possibility for multiple uncomplete.
        all_subdir_temp = []
        pts_to_delete = []
        deleted_points = 0
        
        # sweep all points
        for i,j in enumerate(all_points):
            temp = [k for k in all_subdir if j in k] # all addresses given one point condition
            if len(temp) != num_rep: # compare number of repetition for each point
                                     # if different from the reference (num_rep)
                                     # I will remove this point
                pts_to_delete.append(j)
                
        # removing uncomplete points from the list all_subdir
        if len(pts_to_delete) != 0:
            for i,w in enumerate(all_subdir): # sweep all points
                # ~ print(w)
                for j,x in enumerate(pts_to_delete): # sweep all delete conditions
                    # ~ print(x)
                    if x not in w: # if the current point is not complete, delete from all_subdir
                        all_subdir_temp.append(w)
                    else:
                        deleted_points += 1
                        # print('===== ^ REMOVED ^ =====')
            all_subdir = all_subdir_temp
        del all_subdir_temp
                    
        print('Points deleted because they were not complete',pts_to_delete,'  '+str(deleted_points)+' pt(s)')
        print('Total number of data directories',len(all_subdir))
    else:
        print('No points deleted because they were not complete')
        print('Total number of data directories',len(all_subdir))        

    # all_subdir ce sont tous les répertoires contenant des donnés à analyser
    
    file_cfg = [file_path, dir_path, work_rep, filename]
    slash_cfg = [myslashpos, slashcond]
        
    return file_cfg, slash_cfg, all_subdir
    
def simu_conditions(all_subdir, myslashpos, slashcond,filename):

    # All points of simulation
    all_points = [point[myslashpos[slashcond] + 1:myslashpos[slashcond + 1]] for point in all_subdir]
    all_points = list(dict.fromkeys(all_points))
    print(all_points)

    # Name of the conditions
    condition_separator_position = [m.start() for m in re.finditer('_', all_points[0])]
    if len(condition_separator_position) == 0:
        condition_name = [re.sub('[0-9]+','',all_points[0])]
        condition_numb = [re.sub('[a-zA-Z]','',kk) for _,kk in enumerate(all_points)]
        points_and_coord = [k for j,k in enumerate(all_points)]
        points_and_coord = {y: condition_numb[x] for x,y in enumerate(all_points)}

    else:
        condition_name = [[] for k in range(len(condition_separator_position) + 1)]
        for k, m in enumerate(condition_separator_position):
            condition_name[k] = re.sub('[0-9]+', '', all_points[0][m - condition_separator_position[0]:m])
        condition_name[-1] = re.sub('[0-9]+', '', all_points[0][m + 1:])
        
        # Put together points with their coordinates        
        points_and_coord = {}
        for k,pt in enumerate(all_points):
            print(f'{k:03.0f}','>',pt)
            w = condition_separator_position[0]
            temp = []
            for i,j in enumerate(condition_name):
                temp_cond_pos = pt.find(j)
                try:
                    temp_sep_pos = condition_separator_position[i]
                    temp_cond_num = re.sub(j,'',pt[temp_cond_pos:temp_sep_pos])
                except:
                    temp_cond_num = re.sub(j,'',pt[temp_cond_pos:])
                temp.append(temp_cond_num)
            points_and_coord[pt] = temp
        # {'DC01_RF08': ['01', '08'], 'DC01_RF09': ['01', '09'], 'DC01_RF10': ['01', '10'], ... }

    print('> condition names',condition_name)
    print('> number of points',len(all_points))

    # Conditions de la simulation
    nions_pos = filename.find('_N')            # position nombre d'ions dans filename
    N = int(filename[nions_pos+2:nions_pos+7]) # nombre d'ions
    e_GMol = 50   # energie GMol en eV

    print(f'> N_ions = {N}')
    print(f'> e_GMol = {e_GMol}')

    condition_parameters = [condition_name, N, e_GMol]

    return points_and_coord, condition_parameters

# Recovering two varying condition data
def data_retrieve(all_subdir,points_and_coord, condition_parameters, slash_cfg,**kwargs):
    
    # Récupérer les données pour chaque simu
    # Exécution en quelques minutes
    
    # data0, data2, data4 : nom des fichiers de fluo et température
    # PMvar, Tvar         : variation de fluo et température entre passage Gmol et fin de simu
    # deltaEc, deltaEcRel : variation énergie GMol
    # t_c                 : temps de vol de la GMol entre son apparition (z=-1.5mm) et sa disparition (z=1.5mm)
    # r_LC, v_LC, a_LC    : position, vitesse, accélération des ions Ca+
    # r_LC_clip           : même chose mais en dégageant les ions perdus et beaucoup trop loin (hors pièges)
    # dim_nu              : dimensions x, y et z, nuage Ca+ basé sur r_LC_clip

    myslashpos = slash_cfg[0]
    slashcond = slash_cfg[1]

    # determining number of elements on each repetition
    try:
        num_runs = ['Try'+str(kwargs.get('forcenumtry'))]
        num_runs_aux = [runs[myslashpos[slashcond+1]+1:] for runs in all_subdir if list(points_and_coord.keys())[0] in runs]
        num_runs_aux = list(dict.fromkeys(num_runs_aux))
        num_runs_aux = len(num_runs_aux)
        all_subdir = [x for w,x in enumerate(all_subdir) if (w-kwargs.get('forcenumtry'))%num_runs_aux==0]
    except:
        num_runs = [runs[myslashpos[slashcond+1]+1:] for runs in all_subdir if list(points_and_coord.keys())[0] in runs]
        num_runs = list(dict.fromkeys(num_runs))

    # number of repetitions
    print('> Points |',len(points_and_coord))
    print('> Simulations pour chaque point |', num_runs)
    
    data2 = [[] for i in range(len(points_and_coord))] # file path to SimuType0
    data4 = [[] for i in range(len(points_and_coord))] # file path to SimuType2
    data0 = [[] for i in range(len(points_and_coord))] # file path to SimuType4
    data_address = [[] for i in range(len(points_and_coord))]

    # Variables à deux coordonnées : [point, try]
    shapevar = (len(points_and_coord),len(num_runs))

    PMvar = np.zeros(shapevar)
    Tvar = np.zeros(shapevar)
    deltaEc = np.zeros(shapevar)
    deltaEcRel = np.zeros(shapevar)
    SNR = np.zeros(shapevar)
    t_c = np.zeros(shapevar)
    # r_LC_clip = [[[] for i in range(elem_2)] for j in range(elem_0)]
    dim_nu=zeros((len(points_and_coord),len(num_runs),3))

    t0 = time.clock()
    print("Hello")

    # write variables

    # all files to stat
    # ~ fileload = [[[] for w in range(elem_1)] for i in range(elem_0)]
    
    
    for k, address in enumerate(all_subdir):

    # in-loop variables
        pnt = k // len(num_runs)  # actual point
        rep = k  % len(num_runs)          # actual repetition

        # get only .dat files in each simulation directory
        onlyfiles = [f for f in listdir(address) if isfile(join(address, f)) and not "xva" in f and ".dat" in f]
        # build path file
        data2[pnt].append('{}/{}'.format(address,sort(onlyfiles)[0].strip('.dat')))
        data4[pnt].append('{}/{}'.format(address,sort(onlyfiles)[1].strip('.dat')))
        data0[pnt].append('{}/{}'.format(address,sort(onlyfiles)[2].strip('.dat')))
        data_address[pnt].append(address)

        # load fluorescence and T
        try:
            PMvar[pnt][rep] = find_PM_variation_FinalT(address,sort(onlyfiles)[2].strip('.dat'))[0]
            Tvar[pnt][rep]  = find_PM_variation_FinalT(address,sort(onlyfiles)[2].strip('.dat'))[1]
            SNR[pnt][rep] = find_PM_variation_FinalT(address,sort(onlyfiles)[2].strip('.dat'))[2]
        except:
            PMvar[pnt][rep] = None
            Tvar[pnt][rep]  = None
            SNR[pnt][rep]   = None

        # load Ec variation for GMol
        try:
            deltaEc[pnt][rep] = energy_lost(address,'xva'+sort(onlyfiles)[1][4:].strip('.dat'))[2]
            deltaEcRel[pnt][rep] = energy_lost(address,'xva'+sort(onlyfiles)[1][4:].strip('.dat'))[3]
            t_c[pnt][rep] = energy_lost(address,'xva'+sort(onlyfiles)[1][4:].strip('.dat'))[4]
        except:
            deltaEc[pnt][rep] = None
            deltaEcRel[pnt][rep] = None
            t_c[pnt][rep] = None

        # load cloud size before injection
        try:
            my_file = '{}/xva{}'.format(address,sort(onlyfiles)[2].strip('.dat')[4:])
            r_LC,v_LC,a_LC,v_LC_avg = load_xyz_init_bin_DP(my_file)        

            # filter lost ions
            x_LC_clip = [r_LC[0,x] for x in range(len(r_LC[0,:])) if abs(r_LC[0,x]) <1] #6e-2
            y_LC_clip = [r_LC[1,x] for x in range(len(r_LC[1,:])) if abs(r_LC[1,x]) <1]
            z_LC_clip = [r_LC[2,x] for x in range(len(r_LC[2,:])) if abs(r_LC[2,x]) <1]
    #         r_LC_clip[k][j][:] = [x_LC_clip,y_LC_clip,z_LC_clip]
            r_LC_clip = [x_LC_clip,y_LC_clip,z_LC_clip]
            dim_nu[pnt][rep] = [max(r_LC_clip[l][:]) for l in range(3)]
        except:
            dim_nu[pnt][rep] = None
        
        if not(rep % len(num_runs)):
            print( "Point n°", pnt )
        
        print(f'{pnt:02}','-',f'{rep:02}',' > ',data0[pnt][rep])
        
        if k == len(points_and_coord)*len(num_runs) - 1:
            break
			
    t1 = time.clock() - t0
    print("Time elapsed: ", t1, 's') # CPU seconds elapsed (floating point)
    print("Time elapsed: ", t1/60, 'm') # CPU seconds elapsed (floating point)
    
    data_name = [data_address, data0, data2, data4]
    PMandT = [PMvar, Tvar]
    Gmol_data = [deltaEc, deltaEcRel, SNR, t_c]
    
    return data_name, num_runs, PMandT, Gmol_data, r_LC_clip, dim_nu
    
def data_retrieve_onefile(all_subdir,points_and_coord, condition_parameters, slash_cfg):
    
    # Récupérer les données pour chaque simu
    # Exécution en quelques minutes
    
    # data0, data2, data4 : nom des fichiers de fluo et température
    # PMvar, Tvar         : variation de fluo et température entre passage Gmol et fin de simu
    # deltaEc, deltaEcRel : variation énergie GMol
    # t_c                 : temps de vol de la GMol entre son apparition (z=-1.5mm) et sa disparition (z=1.5mm)
    # r_LC, v_LC, a_LC    : position, vitesse, accélération des ions Ca+
    # r_LC_clip           : même chose mais en dégageant les ions perdus et beaucoup trop loin (hors pièges)
    # dim_nu              : dimensions x, y et z, nuage Ca+ basé sur r_LC_clip

    myslashpos = slash_cfg[0]
    slashcond = slash_cfg[1]

    # determining number of elements on each repetition
    num_runs = [runs[myslashpos[slashcond+1]+1:] for runs in all_subdir if list(points_and_coord.keys())[0] in runs]
    num_runs = list(dict.fromkeys(num_runs))

    # number of repetitions
    print('> Points |',len(points_and_coord))
    print('> Simulations pour chaque point |', num_runs)
    
    data0 = [[] for i in range(len(points_and_coord))] # file path to SimuType0
    data_address = [[] for i in range(len(points_and_coord))]

    # Variables à deux coordonnées : [point, try]
    shapevar = (len(points_and_coord),len(num_runs))

    PMvar = np.zeros(shapevar)
    Tvar = np.zeros(shapevar)
    deltaEc = np.zeros(shapevar)
    deltaEcRel = np.zeros(shapevar)
    SNR = np.zeros(shapevar)
    t_c = np.zeros(shapevar)
    # r_LC_clip = [[[] for i in range(elem_2)] for j in range(elem_0)]
    dim_nu=zeros((len(points_and_coord),len(num_runs),3))

    t0 = time.clock()
    print("Hello")

    # write variables

    # all files to stat
    # ~ fileload = [[[] for w in range(elem_1)] for i in range(elem_0)]
    
    for k, address in enumerate(all_subdir):

    # in-loop variables
        pnt = k // len(num_runs)  # actual point
        rep = k  % len(num_runs)          # actual repetition

        # get only .dat files in each simulation directory
        onlyfiles = [f for f in listdir(address) if isfile(join(address, f)) and not "xva" in f and ".dat" in f]
        # build path file
        data0[pnt].append('{}/{}'.format(address,sort(onlyfiles)[0].strip('.dat')))
        data_address[pnt].append(address)

        # load fluorescence and T
        try:
            PMvar[pnt][rep] = find_PM_variation_FinalT(address,sort(onlyfiles)[0].strip('.dat'))[0]
            Tvar[pnt][rep]  = find_PM_variation_FinalT(address,sort(onlyfiles)[0].strip('.dat'))[1]
            SNR[pnt][rep] = find_PM_variation_FinalT(address,sort(onlyfiles)[0].strip('.dat'))[2]
        except:
            PMvar[pnt][rep] = None
            Tvar[pnt][rep]  = None
            SNR[pnt][rep]   = None

        # load cloud size before injection
        try:
            my_file = '{}/xva{}'.format(address,sort(onlyfiles)[0].strip('.dat')[4:])
            r_LC,v_LC,a_LC,v_LC_avg = load_xyz_init_bin_DP(my_file)        

            # filter lost ions
            x_LC_clip = [r_LC[0,x] for x in range(len(r_LC[0,:])) if abs(r_LC[0,x]) <10e-2] #6e-2
            y_LC_clip = [r_LC[1,x] for x in range(len(r_LC[1,:])) if abs(r_LC[1,x]) <10e-2]
            z_LC_clip = [r_LC[2,x] for x in range(len(r_LC[2,:])) if abs(r_LC[2,x]) <1e-0]
    #         r_LC_clip[k][j][:] = [x_LC_clip,y_LC_clip,z_LC_clip]
            r_LC_clip = [x_LC_clip,y_LC_clip,z_LC_clip]
            dim_nu[pnt][rep] = [max(r_LC_clip[l][:]) for l in range(3)]
        except:
            dim_nu[pnt][rep] = None
        
        if not(rep % len(num_runs)):
            print( "Point n°", pnt )
        
        print(f'{pnt:02}','-',f'{rep:02}',' > ',data0[pnt][rep])

    t1 = time.clock() - t0
    print("Time elapsed: ", t1, 's') # CPU seconds elapsed (floating point)
    print("Time elapsed: ", t1/60, 'm') # CPU seconds elapsed (floating point)
    
    data_name = [data_address, data0]
    PMandT = [PMvar, Tvar]
    Gmol_data = [deltaEc, deltaEcRel, SNR, t_c]
    
    return data_name, num_runs, PMandT, Gmol_data, r_LC_clip, dim_nu

def data_retrieve_RFRelax(all_subdir, points_and_coord, condition_parameters, slash_cfg, mode, **kwargs):
    # Récupérer les données pour chaque simu
    # Exécution en quelques minutes

    # data0, data2, data4 : nom des fichiers de fluo et température
    # PMvar, Tvar         : variation de fluo et température entre passage Gmol et fin de simu
    # deltaEc, deltaEcRel : variation énergie GMol
    # t_c                 : temps de vol de la GMol entre son apparition (z=-1.5mm) et sa disparition (z=1.5mm)
    # r_LC, v_LC, a_LC    : position, vitesse, accélération des ions Ca+
    # r_LC_clip           : même chose mais en dégageant les ions perdus et beaucoup trop loin (hors pièges)
    # dim_nu              : dimensions x, y et z, nuage Ca+ basé sur r_LC_clip

    myslashpos = slash_cfg[0]
    slashcond = slash_cfg[1]

    # determining number of elements on each repetition
    try:
        num_runs = [xx for xx in kwargs.get('forcetry')]
        # num_runs = ['Try'+f'{xx:02d}' for xx in range(kwargs.get('forcecond'))]
    except:
        num_runs = [runs[myslashpos[slashcond + 1] + 1:] for runs in all_subdir if
                    list(points_and_coord.keys())[0] in runs]
        num_runs = list(dict.fromkeys(num_runs))

    # number of repetitions
    print('> Points |', len(points_and_coord))
    print('> Simulations pour chaque point |', num_runs)

    data0 = [[] for i in range(len(points_and_coord))]  # file path to SimuType0 or Q
    data6 = [[] for i in range(len(points_and_coord))]  # file path to SimuType6
    data_address = [[] for i in range(len(points_and_coord))]

    # Variables à deux coordonnées : [point, try]
    shapevar = (len(points_and_coord), len(num_runs))

    PMvar = np.zeros(shapevar)
    Tvar = np.zeros(shapevar)
    deltaEc = np.zeros(shapevar)
    deltaEcRel = np.zeros(shapevar)
    SNR = np.zeros(shapevar)
    t_c = np.zeros(shapevar)
    # r_LC_clip = [[[] for i in range(elem_2)] for j in range(elem_0)]
    dim_nu = zeros((len(points_and_coord), len(num_runs), 3))

    all_subdir_aux = []
    for k,l in enumerate(all_subdir):
        try: # if a certain Try is choosen for initialization
            number = kwargs.get('forcetryinit')[0]
            if l[-2:] == number[-2:]:
                all_subdir_aux.append(l)
        except: # if not choose same try as data to plot
            for _,number in enumerate(num_runs):
                if l[-2:] == number[-2:]:
                    all_subdir_aux.append(l)

    print(all_subdir_aux)
    t0 = time.clock()
    print("Hello")

    # write variables
    # all files to stat
    # ~ fileload = [[[] for w in range(elem_1)] for i in range(elem_0)]

    for k, address in enumerate(all_subdir_aux):

        # in-loop variables
        pnt = k // len(num_runs)  # actual point
        rep = k % len(num_runs)  # actual repetition

        # get only .dat files in each simulation directory
        onlyfiles = [f for f in listdir(address) if isfile(join(address, f)) and not "xva" in f and ".dat" in f]
        #         onlyfiles = [f for f in listdir(fileload[k][j]) if isfile(join(fileload[k][j], f))]
        #         onlyfiles = [i for i in onlyfiles if not "xva" in i]            # vire les xva_...
        #         onlyfiles = [i for i in onlyfiles if ".dat" in i]               # ne garde que les .dat
        # build path file
        # ~ print(sort(onlyfiles))
        # ~ print(pnt)
        data0[pnt].append('{}\{}'.format(address, sort(onlyfiles)[1].strip('.dat')))
        data6[pnt].append('{}\{}'.format(address, sort(onlyfiles)[0].strip('.dat')))
        data_address[pnt].append(address)

        # load Ec variation for GMol
        if mode == 'GMol':
            deltaEc[pnt][rep] = energy_lost(address, 'xva' + sort(onlyfiles)[2][4:].strip('.dat'))[2]
            deltaEcRel[pnt][rep] = energy_lost(address, 'xva' + sort(onlyfiles)[2][4:].strip('.dat'))[3]
            t_c[pnt][rep] = energy_lost(address, 'xva' + sort(onlyfiles)[2][4:].strip('.dat'))[4]
        elif mode == 'stopping_power':
            deltaEc[pnt][rep] = energy_lost(address, 'xva' + sort(onlyfiles)[1][4:].strip('.dat'))[2]
            deltaEcRel[pnt][rep] = energy_lost(address, 'xva' + sort(onlyfiles)[1][4:].strip('.dat'))[3]
            t_c[pnt][rep] = energy_lost(address, 'xva' + sort(onlyfiles)[1][4:].strip('.dat'))[4]
        elif mode == 'RFrelax':
            deltaEc[pnt][rep] = None
            deltaEcRel[pnt][rep] = None
            t_c[pnt][rep] = None

        if not (rep % len(num_runs)):
            print("Point n°", pnt)

        print(f'{pnt:02}', '-', f'{rep:02}', ' > ', data0[pnt][rep])

    # print(my_dico[str(cond_two_name)+str(key2)])
    # print(my_dico[str(cond_one_name)+str(key1)])
    t1 = time.clock() - t0
    print("Time elapsed: ", t1, 's')  # CPU seconds elapsed (floating point)
    print("Time elapsed: ", t1 / 60, 'm')  # CPU seconds elapsed (floating point)

    data_name = [data_address, data0, data6]

    return data_name, num_runs

def load_Temp_init_bin_Lan(str_load, flag_print):      

    fid = open(str_load+'.bin',  'rb')   

    a    = fromfile(fid,  int32, 1)        # Read record start tag   
    aux  = fromfile(fid,  int32, 1)   
    junk = fromfile(fid,  int32, 1)        # Read record stop tag   
    N_ions = aux[0]   

    b    = fromfile(fid,  int32, 1)        # Read record start tag   
    aux  = fromfile(fid,  int32, 1)   
    junk = fromfile(fid,  int32, 1)        # Read record stop tag   
    j_save = aux[0]   

    c    = fromfile(fid,  int32, 1)        # Read record start tag   
    aux  = fromfile(fid,  float64, 1)   
    junk = fromfile(fid,  int32, 1)        # Read record stop tag   
    dt_j_save_next = aux[0]   

    d    = fromfile(fid,  int32, 1)        # Read record start tag   
    aux  = fromfile(fid,  float64, 1)   
    junk = fromfile(fid,  int32, 1)        # Read record stop tag   
    eta = aux[0]   

    e    = fromfile(fid,  int32,   1   )        # Read record start tag   
    Temp = fromfile(fid,  float64, 1)   
    junk = fromfile(fid,  int32,   1   )        # Read record stop tag   

    f    = fromfile(fid,  int32,   1   )        # Read record start tag 
    save_T = fromfile(fid,  float64, 7*j_save)   
    junk = fromfile(fid,  int32,   1   )        # Read record stop tag

    fid.close   
    
    print('len save_T', junk)   

    if flag_print == 1:
        print(a, b, c, d, e, f)   
        print('N_ions', N_ions)   
        print('j_save', j_save)   
        print('dt_j_save_next', dt_j_save_next)   
        print('eta', eta)   
        print('Temp', Temp)   


    save_T = reshape(save_T, (j_save, 7), order='F')   
    print(save_T[0])   

    return N_ions, j_save, dt_j_save_next, eta, Temp, save_T 

def data_check(all_subdir,points_and_coord, condition_parameters, slash_cfg):
    
    # Récupérer les données pour chaque simu
    # Exécution en quelques minutes
    
    # data0, data2, data4 : nom des fichiers de fluo et température
    # PMvar, Tvar         : variation de fluo et température entre passage Gmol et fin de simu
    # deltaEc, deltaEcRel : variation énergie GMol
    # t_c                 : temps de vol de la GMol entre son apparition (z=-1.5mm) et sa disparition (z=1.5mm)
    # r_LC, v_LC, a_LC    : position, vitesse, accélération des ions Ca+
    # r_LC_clip           : même chose mais en dégageant les ions perdus et beaucoup trop loin (hors pièges)
    # dim_nu              : dimensions x, y et z, nuage Ca+ basé sur r_LC_clip

    myslashpos = slash_cfg[0]
    slashcond = slash_cfg[1]

    # determining number of elements on each repetition
    num_runs = [runs[myslashpos[slashcond+1]+1:] for runs in all_subdir if list(points_and_coord.keys())[0] in runs]
    num_runs = list(dict.fromkeys(num_runs))

    # number of repetitions
    print('> Points |',len(points_and_coord))
    print('> Simulations pour chaque point |', num_runs)
    
    data2 = [[] for i in range(len(points_and_coord))] # file path to SimuType0
    data4 = [[] for i in range(len(points_and_coord))] # file path to SimuType2
    data0 = [[] for i in range(len(points_and_coord))] # file path to SimuType4
    data_address = [[] for i in range(len(points_and_coord))]
    
    for k, address in enumerate(all_subdir):

    # in-loop variables
        pnt = k // len(num_runs)  # actual point
        rep = k  % len(num_runs)          # actual repetition

        # get only .dat files in each simulation directory
        onlyfiles = [f for f in listdir(address) if isfile(join(address, f))]
        
        if len(onlyfiles) != 15:
            print(address)
        
    return

#############################
### Single case functions ###

def file_name_retrieve(address,simu_type):
        # get only .dat files in each simulation directory
    onlyfiles = [f for f in listdir(address) if isfile(join(address, f)) and not "xva" in f and ".dat" in f]
    tmp_file = []
    # build path file
    if simu_type == 'Collision':
        for k in range(3):
            tmp_file.append('{}/{}'.format(address,sort(onlyfiles)[k].strip('.dat')))
    elif simu_type == 'RF_relax':
        for k in range(4):
            try:
                tmp_file.append('{}/{}'.format(address,sort(onlyfiles)[k].strip('.dat')))
            except(IndexError):
                tmp_file.append('No file for this SimuType')
                print('A file is missing, skiping the file number',k)
            
    return onlyfiles, tmp_file

# Get the size of the cloud
def cloud_size(address,onlyfiles, **kwargs):
    radial_thresh = kwargs.get('radial_thresh', 10e-2)
    axial_thresh  = kwargs.get('axial_thresh', 1e-0)
    
    try:
        my_file = '{}/xva{}'.format(address,sort(onlyfiles)[0].strip('.dat')[4:])
        r_LC,v_LC,a_LC = load_xyz_init_bin_DP(my_file)

        # filter lost ions
        x_LC_clip = [r_LC[0,x] for x in range(len(r_LC[0,:])) if abs(r_LC[0,x]) < radial_thresh] #6e-2
        y_LC_clip = [r_LC[1,x] for x in range(len(r_LC[1,:])) if abs(r_LC[1,x]) < radial_thresh]
        z_LC_clip = [r_LC[2,x] for x in range(len(r_LC[2,:])) if abs(r_LC[2,x]) < axial_thresh]
#         r_LC_clip[k][j][:] = [x_LC_clip,y_LC_clip,z_LC_clip]
        r_LC_clip = [x_LC_clip,y_LC_clip,z_LC_clip]
        dim_nu = [max(r_LC_clip[l][:]) for l in range(3)]
    except:
        r_LC_clip = None
        dim_nu = None
            
    return r_LC_clip, dim_nu

# Return Temp, Fluo, for SimuType2 RF with laser
def load_T_PM_cloud_GMol(address,flag_plot):
    # retrieve files adresses
    onlyfiles, data_files = file_name_retrieve(address,simu_type='Collision')
    file0 = data_files[0];file2 = data_files[1];file4 = data_files[2]
    # load Ec variation for GMol
    _,_,deltaEc,deltaEcRel,t_c = energy_lost(address,'xva'+sort(onlyfiles)[2][4:].strip('.dat'))
    # load cloud size before injection
    r_LC_clip, dim_nu = cloud_size(address,onlyfiles)
    # load fluorescence and T var
    # Temperature during time
    temps, temperature, fluo = plot_T_and_PM_Init_Inje_Evol(address+'/',sort(onlyfiles)[0].strip('.dat')[4:],flag_plot,'T_PM',xlim1=(-1.5,55),ylim1=(5e-5,12e3))

    fluo_var = [fluo[1], temperature[1], fluo[2]]                        # temps = [tt,t_aux1, t_aux2]
    GMol_var = [deltaEc,deltaEcRel,t_c]                                  # fluo = [PM, PM_variation, SNR]
    cloud_atlas = [temps[0], temperature[0], temperature[2], fluo[0]]    # temperature = [T_aux, T_variation]
    t_aux = [temps[1],temps[2]]
    
    return fluo_var, GMol_var, cloud_atlas, t_aux, r_LC_clip

# Return Temp, Fluo, for SimuType6 RF without laser
def load_RF_relax(address,flag_plot,RF_relax_type):
    # retrieve files adresses
    onlyfiles, data_files = file_name_retrieve(address,simu_type='RF_relax')
    file0 = data_files[0];file2 = data_files[1];file4 = data_files[2];file6 = data_files[3]

    # load fluorescence and T var
    # Temperature during time
    if RF_relax_type == 'AfterCool':
        tt, T_aux, T_CM = plot_T_and_PM_InitQ_Evol_AfterCool(address+'/',flag_plot,'RFrelax_AfterCool',xlim1=(-1.5,55),ylim1=(5e-5,12e3))
        GMol_var = None
        r_LC_clip = None
    elif RF_relax_type == 'AfterInj':
        # load Ec variation for GMol
        _, _, deltaEc, deltaEcRel, t_c = energy_lost(address, 'xva' + sort(onlyfiles)[2][4:].strip('.dat'))
        # load cloud size before injection
        r_LC_clip, dim_nu = cloud_size(address, onlyfiles)
        tt, T_aux, T_CM = plot_T_and_PM_Init_RFrelax_AfterInj_Evol(address+'/',sort(onlyfiles)[0].strip('.dat')[4:],flag_plot,'RFrelax_AfterInj',xlim1=(-1.5,55),ylim1=(5e-5,12e3))
        GMol_var = [deltaEc, deltaEcRel, t_c]

    return tt, T_aux, T_CM, GMol_var, r_LC_clip

# ======================================================================== #
# pour les couleurs
from matplotlib import colors as mcolors
def cc(arg,alpha):
    '''
    Shorthand to convert 'named' colors to rgba format at x% opacity.
    '''
    return mcolors.to_rgba(arg, alpha=alpha)

# convert H:M:s time to sec
def get_sec(time_str):
    h, m, s = time_str.split(':')
#     return h,m,s
    return int(h) * 3600 + int(m) * 60 + int(s)
