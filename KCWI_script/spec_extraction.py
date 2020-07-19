# !/usr/bin/python3
# -*- coding: utf-8 -*-
# @Time      :   2020/7/8 2:04 PM
# @Author    :   Shiwu
# @Site      :   
# @File      :   spec_extraction.py
# @Desc      :   
# @Software  :   PyCharm
# @license   :   Copyright(C), C
# @Contact   :   shiwuzhang95@gmail.com


'''
extract the spectra from some windows
'''

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const
from scipy.interpolate import interp2d
from Class_func import Map

from Map import lyemissioncube,heiiemissioncube,civemissioncube ,ly_noise,heii_noise,civ_noise,lyvelocitymap,heiivelocitymap,civvelocitymap

class spec_extractor:

    def __init__(self,cube=None):
        self.cube=cube
    def snr_img(self,cube,noise,threshold,scale=1e-17):

        nb_img=cube.optimal_img*scale
        noise_stacked=noise*scale*np.sqrt(cube.shape[0])
        snr=nb_img/noise_stacked
        snr[snr<=threshold]=threshold

        return snr

    def img_trim(self,img,trim_edge,fill_value):

        trim_bool=np.ones_like(img).astype(np.bool)
        trim_bool[trim_edge[0]:trim_edge[1],trim_edge[2]:]=False
        img[trim_bool]=fill_value

        return  img

    def img_interpolate(self,x,y,xnew,ynew,img):

        interfunc=interp2d(y,x,img,kind='linear')
        img_new=interfunc(ynew,xnew)

        return img_new

    def cube_interpolate(self,x,y,xnew,ynew,cube):

        cube_new=np.zeros((cube.shape[0],xnew.size,ynew.size))
        for i in range(cube.shape[0]):
            cube_new[i,:,:]=self.img_interpolate(x,y,xnew,ynew,cube[i,:,:])

        return cube_new

    def img_redirection(self,img):

        img=np.rot90(img)[:,::-1]

        return img

    def cube_redirection(self,cube):

        for i in range(cube.shape[0]):
            cube[i,:,:]=self.img_redirection(cube[i,:,:])

        return cube

    def radius_mesh(self,x,y):

        x_mesh,y_mesh=np.meshgrid(y,x)
        r_mesh=np.sqrt((x_mesh**2)+(y_mesh**2))

        return r_mesh

    def spec_extract(self,cube,snr1,snr2,snr3,r_mesh,
                     nx,ny,delta_x,delta_y,thre1,thre2):

        lib,x_lib,y_lib,r_lib=[],[],[],[]

        for i in range(nx, cube.shape[1] - nx, nx + delta_x):
            for j in range(ny, cube.shape[2] - ny, ny + delta_y):
                cell_cube = cube[:, i - nx:i + nx, j - ny:j + ny]
                cell_snr1 = snr1[i - nx:i + nx, j - ny:j + ny]
                cell_snr2 = snr2[i - nx:i + nx, j - ny:j + ny]
                cell_snr3 = snr3[i - nx:i + nx, j - ny:j + ny]
                cell_cube1 = cell_cube[:, (cell_snr1 > thre1) & (cell_snr2 > thre1) & (cell_snr3 > thre2)]

                if cell_cube1.size > 0:
                    spec = cell_cube.mean(axis=1).mean(axis=1)
                    r = r_mesh[i - nx:i + nx, j - ny:j + ny]
                    lib.append(spec)
                    x_lib.append(i)
                    y_lib.append(j)
                    r_lib.append(r.min())

        return np.array(lib), np.array(x_lib), np.array(y_lib), np.array(r_lib)


def snr(cube_lib,noise_lib):

    extractor=spec_extractor()

    thre_lib=[4.,3.,3.]
    snr_lib=[]
    for i in range(len(cube_lib)):
        snr_img=extractor.snr_img(cube_lib[i],noise_lib[i],thre_lib[i])
        snr_lib.append(snr_img)

    return snr_lib

def interpolate(img_lib,cube_lib,v_lib,x,y):

    xnew=np.linspace(x.min(),x.max(),400)
    ynew = np.linspace(y.min(), y.max(), 400)

    extractor=spec_extractor()
    for i in range(len(img_lib)):
        img_lib[i]=extractor.img_interpolate(x,y,xnew,ynew,img_lib[i])
        cube_lib[i]=extractor.cube_interpolate(x,y,xnew,ynew,cube_lib[i])
        v_lib[i]=extractor.img_interpolate(x,y,xnew,ynew,v_lib[i])

    return img_lib,cube_lib,v_lib,xnew,ynew

def redirection(img_lib,cube_lib,v_lib):

    extractor=spec_extractor()
    for i in range(len(img_lib)):
        img_lib[i]=extractor.img_redirection(img_lib[i])
        cube_lib[i] = extractor.cube_redirection(cube_lib[i])
        v_lib[i]=extractor.img_redirection(v_lib[i])

    return img_lib,cube_lib,v_lib

ly_dataedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':4008*u.AA,'zhi':4039*u.AA}#4002,4040
ly_noiseedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':4050*u.AA,'zhi':4100*u.AA}
ly_continuumedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':4100*u.AA,'zhi':4300*u.AA}
lymap=Map(filename='/Users/shiwuzhang/WS/ASTRO/MAMMOTH_KCWI/1441+4003_comb_ss_icubes.fits',
          refpoint=np.array([220.3520886,40.05269183])*u.deg,dataedge=ly_dataedge,
          continuumedge=ly_continuumedge,noiseedge=ly_noiseedge,
          noiselevel=2)
_,_,_,lyemissioncube2=lymap.optimalmap()

heii_dataedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5455*u.AA,'zhi':5496*u.AA}#5465,5485
heii_noiseedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5400*u.AA,'zhi':5450*u.AA}
heii_continuumedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5500*u.AA,'zhi':5600*u.AA}
heiimap=Map(filename='/Users/shiwuzhang/WS/ASTRO/MAMMOTH_KCWI/1441+4003_comb_psfs_icubes.fits',
            refpoint=np.array([220.3520886, 40.05269183]) * u.deg,dataedge=heii_dataedge,
            continuumedge=heii_continuumedge,noiseedge=heii_noiseedge,
            noiselevel=1.2)#1
_,_,_,heiiemissioncube2=heiimap.optimalmap()

civ_dataedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5154*u.AA,'zhi':5198*u.AA}#5203
civ_noiseedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5100*u.AA,'zhi':5150*u.AA}
civ_continuumedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5000*u.AA,'zhi':5100*u.AA}
civmap=Map(filename='/Users/shiwuzhang/WS/ASTRO/MAMMOTH_KCWI/1441+4003_comb_psfs_icubes.fits',
           refpoint=np.array([220.3520886, 40.05269183]) * u.deg,dataedge=civ_dataedge,
           continuumedge=civ_continuumedge,noiseedge=civ_continuumedge,
           noiselevel=1.2)#1.5
_,_,_,civemissioncube2=civmap.optimalmap()


extractor=spec_extractor()

emissioncube=[lyemissioncube,heiiemissioncube,civemissioncube]
noise=[ly_noise,heii_noise,civ_noise]
snr_lib=snr(emissioncube,noise)
snr_lib[1]=extractor.img_trim(snr_lib[1],[25,44,9],3.)
snr_lib[2]=extractor.img_trim(snr_lib[2],[23,46,8],3.)

lyv=lyvelocitymap.copy()
lyv[np.isnan(lyv)]=0.
heiiv=heiivelocitymap.copy()
heiiv[np.isnan(heiiv)]=0.
civv=civvelocitymap.copy()
civv[np.isnan(civv)]=0.
v_lib=[lyv,heiiv,civv]

x=heiiemissioncube.ra.to(u.arcsec).value
y=heiiemissioncube.dec.to(u.arcsec).value

emissioncube2=[lyemissioncube2._data,heiiemissioncube2._data,civemissioncube2._data]
snr_lib,emissioncube2,v_lib,xnew,ynew=interpolate(snr_lib,emissioncube2,v_lib,x,y)
snr_lib,emissioncube2,v_lib=redirection(snr_lib,emissioncube2,v_lib)


r_mesh_new=extractor.radius_mesh(xnew,ynew)
result_spec_ly=extractor.spec_extract(emissioncube2[0],snr_lib[1],
                                      snr_lib[2],snr_lib[0],r_mesh_new,10,10,13,13,3.5,12.)
spec_ly_lib,x_window_ly,y_window_ly,r_ly=result_spec_ly

result_spec_ly2=extractor.spec_extract(emissioncube2[0],snr_lib[1],
                                      snr_lib[2],snr_lib[0],r_mesh_new,10,10,13,13,1.,14.)
spec_ly_lib2,x_window_ly2,y_window_ly2,r_ly2=result_spec_ly2

result_spec_heii=extractor.spec_extract(emissioncube2[1],snr_lib[1],
                                        snr_lib[2],snr_lib[0],r_mesh_new,10,10,13,13,3.5,1.)
spec_heii_lib,x_window_heii,y_window_heii,r_heii=result_spec_heii

result_spec_civ=extractor.spec_extract(emissioncube2[2],snr_lib[1],
                                       snr_lib[2],snr_lib[0],r_mesh_new,10,10,13,13,3.5,1.)
spec_civ_lib,x_window_civ,y_window_civ,r_civ=result_spec_civ


lamda_ly=lyemissioncube2.spectral_axis.to(u.AA).value
lamda_heii=heiiemissioncube2.spectral_axis.to(u.AA).value
lamda_civ=civemissioncube2.spectral_axis.to(u.AA).value
lamda0_ly=1215.67*3.3076
lamda0_heii=1640*3.3375
lamda0_civ=1549*3.339

v_ly=((lamda_ly - lamda0_ly) / lamda0_ly) * const.c.to(u.km / u.s).value
v_heii=((lamda_heii - lamda0_heii) / lamda0_heii) * const.c.to(u.km / u.s).value
v_civ=((lamda_civ - lamda0_civ) / lamda0_civ) * const.c.to(u.km / u.s).value