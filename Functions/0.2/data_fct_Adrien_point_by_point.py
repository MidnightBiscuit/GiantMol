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
    v_LC = 1.e3*aux[3:6,0:n_ions1]

    #xyz = 1.e3*aux[0:3,:]
    a_LC = 1.e3*aux[6:9,0:n_ions1]

    return r_LC,v_LC,a_LC

def plot_XYZ(file_name,fig_name='2',fig_title='XYZ'):
    r_LC,v_LC,a_LC = load_xyz_init_bin_DP(file_name)
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
    file1 = 'SimuType0'    + file_name[i_aux:] ########## file1 = 'SimuType0'    + file_name[i_aux:17+36]
    file2 = 'SimuType4_01' + file_name[i_aux:]
    file3 = 'SimuType2_01' + file_name[i_aux:]

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
    T_aux = concatenate( (T_aux1,T_aux2,T_aux3) )
    PM    = concatenate( (PM1,PM2,PM3) )
    
    if flag_plot == 1 :
        #fig_name = file_name[-9:]
        figure(fig_name); clf()
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
        plt.setp(ax1.get_xticklabels(), visible=False)

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
        
    return tt, T_aux, PM, PM_variation, T_variation

def find_PM_variation_FinalT(file_dir2,file_name):
    i_aux = file_name.find('_N')+1
    file1 = 'SimuType0_{}'.format(file_name[i_aux:].strip('50eV.dat'))    # for N=1024
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
        if 'Hobitton' in file_path :
            slashcond = 5 # the slash number just before runs (after date)
        else:
            slashcond = 6 # the slash number just before runs (after date)
    else:
        slashcond = 2

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
        for i,j in enumerate(all_points):
            temp = [k for k in all_subdir if j in k]
            if len(temp) != num_rep: # compare number of repetition for each point
                                     # if different from the reference (num_rep)
                                     # I will remove this point
                pts_to_delete.append(j)
        # removing uncomplete points from the list all_subdir
        for i,w in enumerate(all_subdir): # sweep all points
            for j,x in enumerate(pts_to_delete): # sweep all delete conditions
				# print(w)
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

    # Name of the conditions
    condition_separator_position = [m.start() for m in re.finditer('_', all_points[0])]
    condition_name = [[] for k in range(len(condition_separator_position) + 1)]

    for k, m in enumerate(condition_separator_position):
        condition_name[k] = re.sub('[0-9]+', '', all_points[0][m - condition_separator_position[0]:m])
    condition_name[-1] = re.sub('[0-9]+', '', all_points[0][m + 1:])
    print('> condition names',condition_name)

    print('> number of points',len(all_points))

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

    # Conditions de la simulation
    nions_pos = filename.find('_N')            # position nombre d'ions dans filename
    N = int(filename[nions_pos+2:nions_pos+7]) # nombre d'ions
    e_GMol = 50   # energie GMol en eV

    print(f'> N_ions = {N}')
    print(f'> e_GMol = {e_GMol}')

    condition_parameters = [condition_name, N, e_GMol]

    return points_and_coord, condition_parameters

# Recovering two varying condition data

def data_retrieve(all_subdir,points_and_coord, condition_parameters, slash_cfg):
    
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
    data2 = [[] for i in range(len(points_and_coord))] # file path to SimuType2
    data4 = [[] for i in range(len(points_and_coord))] # file path to SimuType4
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

    # ~ fileload = [[[] for w in range(elem_1)] for i in range(elem_0)]
    # ~ for rf in range(elem_2):     # Vrf  j
        # ~ for dc in range(elem_1):
            # ~ print('rf - dc - tr')
            # ~ for tr in range(elem_0): # try  k
                # ~ address = ( work_rep[:myslashpos[slashcond]]
                            # ~ +str(cond_zero_name)+str(tr)
                            # ~ +str(cond_one_name)+str(key1[dc])
                            # ~ +str(cond_two_name)+str(key2[rf]) )
                # ~ fileload[tr][dc].append(adress)
                # ~ print(f'{rf:02}','-',f'{dc:02}','-',f'{tr:02}',' > ',fileload[tr][dc][rf])

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
#         onlyfiles = [f for f in listdir(fileload[k][j]) if isfile(join(fileload[k][j], f))]
#         onlyfiles = [i for i in onlyfiles if not "xva" in i]            # vire les xva_...
#         onlyfiles = [i for i in onlyfiles if ".dat" in i]               # ne garde que les .dat
        # build path file
        data0[pnt].append('{}/{}'.format(address,sort(onlyfiles)[0].strip('.dat')))
        data2[pnt].append('{}/{}'.format(address,sort(onlyfiles)[1].strip('.dat')))
        data4[pnt].append('{}/{}'.format(address,sort(onlyfiles)[2].strip('.dat')))
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

        # load Ec variation for GMol
        deltaEc[pnt][rep] = energy_lost(address,'xva'+sort(onlyfiles)[2][4:].strip('.dat'))[2]
        deltaEcRel[pnt][rep] = energy_lost(address,'xva'+sort(onlyfiles)[2][4:].strip('.dat'))[3]
        t_c[pnt][rep] = energy_lost(address,'xva'+sort(onlyfiles)[2][4:].strip('.dat'))[4]

        # load cloud size before injection
        try:
            my_file = '{}/xva{}'.format(address,sort(onlyfiles)[0].strip('.dat')[4:])
            r_LC,v_LC,a_LC = load_xyz_init_bin_DP(my_file)        

            # filter lost ions
            x_LC_clip = [r_LC[0,x] for x in range(len(r_LC[0,:])) if abs(r_LC[0,x]) <6e-2]
            y_LC_clip = [r_LC[1,x] for x in range(len(r_LC[1,:])) if abs(r_LC[1,x]) <6e-2]
            z_LC_clip = [r_LC[2,x] for x in range(len(r_LC[2,:])) if abs(r_LC[2,x]) <1e-0]
    #         r_LC_clip[k][j][:] = [x_LC_clip,y_LC_clip,z_LC_clip]
            r_LC_clip = [x_LC_clip,y_LC_clip,z_LC_clip]
            dim_nu[pnt][rep] = [max(r_LC_clip[l][:]) for l in range(3)]
        except:
            dim_nu[pnt][rep] = None
        
        if not(rep % len(num_runs)):
            print( "Point n°", pnt )
        
        print(f'{pnt:02}','-',f'{rep:02}',' > ',data0[pnt][rep])

    # print(my_dico[str(cond_two_name)+str(key2)])
    # print(my_dico[str(cond_one_name)+str(key1)])
    t1 = time.clock() - t0
    print("Time elapsed: ", t1, 's') # CPU seconds elapsed (floating point)
    print("Time elapsed: ", t1/60, 'm') # CPU seconds elapsed (floating point)
    
    data_name = [data_address, data0, data2, data4]
    PMandT = [PMvar, Tvar]
    Gmol_data = [deltaEc, deltaEcRel, SNR, t_c]
    
    return data_name, num_runs, PMandT, Gmol_data, r_LC_clip, dim_nu

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
