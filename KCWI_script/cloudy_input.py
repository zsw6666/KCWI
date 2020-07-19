# !/usr/bin/python3
# -*- coding: utf-8 -*-
# @Time      :   2020/7/3 9:21 AM
# @Author    :   Shiwu
# @Site      :   
# @File      :   run_cloudy.py
# @Desc      :   
# @Software  :   PyCharm
# @license   :   Copyright(C), C
# @Contact   :   shiwuzhang95@gmail.com


import numpy as np
import pandas as pd
import pyCloudy as pc
import os


pc.config.cloudy_exe='/home/zsw666/cloudy/c17.01/source/cloudy.exe'

def make_model(dic_,model_name,dens,ab0,n_column,ionpar):

    full_model_name='{0}_nH{1:.2f}_Z{2:.1f}_NH{3:.1f}_U{4:.2f}'.format(model_name,dens,ab0,n_column,ionpar)
    emis_tab=['H  1     1215.67A','h  1     6562.81A','c  4  1548.19A','he 2     1640.43A']
    c_input=pc.CloudyInput('{0}{1}'.format(dic_,full_model_name))
    c_input.set_star(SED='table agn',SED_params='6.0 -1.4 -0.5 -1.0',lumi_unit='ionization parameter',lumi_value=ionpar)
    c_input.set_sphere(sphere=False)
    c_input.set_radius(r_in=23,r_out=23.5)
    c_input.set_distance(dist=18544.,unit='Mpc',linear=True)
    c_input.set_abund(predef='GASS',nograins=True,metals=ab0)
    c_input.set_cste_density(dens)
    c_input.set_other(('covering factor 0.3','CMB redshift 2.3116','set trim -20',
                       'turbulence 50 km/s','COSMIC RAY BACKGROUND','print last'))
    c_input.set_emis_tab(emis_tab)
    c_input.set_iterate(to_convergence=True)
    c_input.set_stop(stop_criter=['temperature 10K linear','column density %d'%n_column,'neutral column density 17.2'])
    c_input.print_input(to_file=True,verbose=False)


def correct_input(dic,modelname):
    file_list = [file for file in os.listdir(dic) if modelname in file]

    for file in file_list:

        with open(dic + file, 'r+') as f:
            new_f = f.readlines()
            f.seek(0)
            for line in new_f:
                if 'filling factor' not in line:
                    f.write(line)
            f.truncate()
dic='/home/zsw666/cloudy/model_photo/'
model_name='CGM'
tab_density=np.array([0.1,1,10])
tab_abs=np.linspace(0.1,5,10)
N_column=np.linspace(18,22,10)
ionparameters=np.linspace(-3,0,5)
for dens in tab_density:
    for abund in tab_abs:
        for nc in N_column:
            for u in ionparameters:
                make_model(dic_=dic,model_name=model_name,dens=dens,ab0=abund,n_column=nc,ionpar=u)

correct_input(dic,model_name)

#pc.run_cloudy(dir_=dic,model_name=model_name,use_make=True,n_proc=20)
