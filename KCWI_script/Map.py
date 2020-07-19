# !/usr/bin/python3
# -*- coding: utf-8 -*-
# @Time      :   2019/12/10 9:45 PM
# @Author    :   Shiwu
# @Site      :   
# @File      :   Map.py
# @Desc      :   
# @Software  :   PyCharm
# @license   :   Copyright(C), C
# @Contact   :   shiwuzhang95@gmail.com

'''
This script is to generate flux map, velocity map and
dispersion map of KCWI data, this is the MAMMOTH-1 field
'''



'''
这个脚本读取fits，利用cube得到optimal-extraction和
velocity map, dispersion map。
主要步骤是：
    1.  从原始cube产生3种emission, noise, continuum这三种cube，
    emission cube截取发射线部分，noise cube和continuum cube避开发射线，
    因为这两种cube分别是用来计算noise和continuum的。
    2.对continuum cube做median stack产生continuum image，之后做
    emission-continuum for each slice of emission cube。
    3.for each slice of emission cube,认为>n*sigma的部分是信号，
    n是noise level, sigma是standard deviation（从noise cube中计算得到），mask
    掉不是信号的部分产生mask cube --> mask cube*emission cube=optimal cube.
    4. 计算velocity map dispersion map and SNR map.
'''


import numpy as np
from astropy import units as u
from Class_func import Map, Img_interpsmooth
import matplotlib.pyplot as plt


#定义每一种cube的边界
ly_dataedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':4002*u.AA,'zhi':4040*u.AA}#4002,4040
ly_noiseedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':4050*u.AA,'zhi':4100*u.AA}
ly_continuumedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':4100*u.AA,'zhi':4300*u.AA}
heii_dataedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5465*u.AA,'zhi':5485*u.AA}#5465,5485
heii_noiseedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5400*u.AA,'zhi':5450*u.AA}
heii_continuumedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5500*u.AA,'zhi':5600*u.AA}
civ_dataedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5163*u.AA,'zhi':5183*u.AA}#5203
civ_noiseedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5100*u.AA,'zhi':5150*u.AA}
civ_continuumedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5000*u.AA,'zhi':5100*u.AA}

lymap=Map(filename='/Users/shiwuzhang/WS/ASTRO/MAMMOTH_KCWI/1441+4003_comb_ss_icubes.fits',
          refpoint=np.array([220.3520886,40.05269183])*u.deg,dataedge=ly_dataedge,
          continuumedge=ly_continuumedge,noiseedge=ly_noiseedge,
          noiselevel=2)

heiimap=Map(filename='/Users/shiwuzhang/WS/ASTRO/MAMMOTH_KCWI/1441+4003_comb_psfs_icubes.fits',
            refpoint=np.array([220.3520886, 40.05269183]) * u.deg,dataedge=heii_dataedge,
            continuumedge=heii_continuumedge,noiseedge=heii_noiseedge,
            noiselevel=1.2)#1

civmap=Map(filename='/Users/shiwuzhang/WS/ASTRO/MAMMOTH_KCWI/1441+4003_comb_psfs_icubes.fits',
           refpoint=np.array([220.3520886, 40.05269183]) * u.deg,dataedge=civ_dataedge,
           continuumedge=civ_continuumedge,noiseedge=civ_continuumedge,
           noiselevel=1.2)#1.5

dic='/Users/shiwuzhang/WS/ASTRO/MAMMOTH_KCWI/'
lyoptimalcube,lynoisecube,maskcube,lyemissioncube=lymap.optimalmap()
lysnrmap=lymap.snrmap()
lysnrmap[lysnrmap<2]=0
lysnrmap=lymap.img_cut(0.99,[2,-1,2,-1],lysnrmap)
lyvelocitymap=lymap.momentmap(order=1,redshift= 2.3076,restvalue=1215.67*u.AA)
lydispersionmap=lymap.momentmap(order=2,redshift= 2.3076,restvalue=1215.67*u.AA)
delta=lyoptimalcube.spectral_axis.to(u.AA).value[1]-lyoptimalcube.spectral_axis.to(u.AA).value[0]
lyoptimalcube.optimal_img=lyoptimalcube.optimal_img*10*delta
lyemissioncube.optimal_img=lyemissioncube._data.sum(axis=0)*10*delta
ly_noise=lynoisecube.noise().mean()*10*delta
lydispersionmap[lydispersionmap>1000]=np.nan
lydispersionmap[lydispersionmap<200]=np.nan
lydispersionmap[lysnrmap<2]=np.nan
lyvelocitymap[np.isnan(lydispersionmap)]=np.nan
lyoptimalcube.optimal_img[lysnrmap<2]=0


heiioptimalcube,heiinoisecube,_,heiiemissioncube=heiimap.optimalmap()
heiisnrmap=heiimap.snrmap()
heii_noise=heiinoisecube.noise().mean()*10*delta
heiisnrmap[heiisnrmap<2]=0
heiisnrmap[23:29,9:11]=0
heiisnrmap=heiimap.img_cut(0.95,[18,48,10,-2],heiisnrmap)
heiivelocitymap=heiimap.momentmap(order=1,redshift=2.3375,restvalue=1640*u.AA)
heiidispersionmap=heiimap.momentmap(order=2,redshift=2.3375,restvalue=1640*u.AA)
delta=heiioptimalcube.spectral_axis.to(u.AA).value[1]-heiioptimalcube.spectral_axis.to(u.AA).value[0]
heiioptimalcube.optimal_img=heiioptimalcube.optimal_img*10*delta
heiispec=heiiemissioncube.with_spectral_unit(u.km/u.s,velocity_convention='relativistic',
                                             rest_value=(1+2.3375)*1640*u.AA).spectral_axis.value
heiiemissioncube.optimal_img=heiiemissioncube._data.sum(axis=0)*10*delta
heiidispersionmap[heiidispersionmap>1000]=np.nan
heiidispersionmap[heiidispersionmap<200]=np.nan
heiidispersionmap[heiisnrmap<2]=np.nan
# heiivelocitymap[np.isnan(heiidispersionmap)]=np.nan
heiioptimalcube.optimal_img[heiisnrmap<2]=0


civoptimalcube,civnoisecube,_,civemissioncube=civmap.optimalmap()
civsnrmap=civmap.snrmap()
civ_noise=civnoisecube.noise().mean()*10*delta
civsnrmap[civsnrmap<2]=0
civsnrmap[15:21,7:13]=0
civsnrmap=civmap.img_cut(0.95,[18,50,8,-2],civsnrmap)
civvelocitymap=civmap.momentmap(order=1,redshift=2.339,restvalue=1549*u.AA)#34385
civdispersionmap=civmap.momentmap(order=2,redshift=2.339,restvalue=1549*u.AA)
civoptimalcube.optimal_img=civoptimalcube.optimal_img*10*delta
civspec=civemissioncube.with_spectral_unit(u.km/u.s,velocity_convention='relativistic',
                                           rest_value=(1+2.339)*1549*u.AA).spectral_axis.value
civemissioncube.optimal_img=civemissioncube._data.sum(axis=0)*10*delta
civdispersionmap[civdispersionmap>1000]=np.nan
civdispersionmap[civdispersionmap<200]=np.nan
civdispersionmap[civsnrmap<2]=np.nan
# civvelocitymap[np.isnan(civdispersionmap)]=np.nan
civoptimalcube.optimal_img[civsnrmap<2]=0