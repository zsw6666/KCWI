# !/usr/bin/python3
# -*- coding: utf-8 -*-
# @Time      :   2020/1/19 11:24 AM
# @Author    :   Shiwu
# @Site      :   
# @File      :   test.py
# @Desc      :   
# @Software  :   PyCharm
# @license   :   Copyright(C), C
# @Contact   :   shiwuzhang95@gmail.com


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
from Class_func import Map,Img_interpsmooth
import matplotlib.pyplot as plt


def preprocess(optimalcube, velocitymap, dispersionmap, snrmap, noise_level=3):
    '''
    before plot the image, do
    preprocess to image, optimalcube
    should be an optimalcube with optimal
    image
    '''
    velocitymap[np.isnan(velocitymap)] = 0
    dispersionmap[np.isnan(dispersionmap)] = 0

    optimalimg_inter, x_inter, y_inter = Img_interpsmooth(optimalcube.optimal_img,
                                                          optimalcube.ra.to(u.arcsec).value,
                                                          optimalcube.dec.to(u.arcsec).value, [5, 9])
    snrmap_inter, _, _ = Img_interpsmooth(snrmap,
                                          optimalcube.ra.to(u.arcsec).value,
                                          optimalcube.dec.to(u.arcsec).value, [5, 9])

    velocitymap_inter, _, _ = Img_interpsmooth(velocitymap,
                                               optimalcube.ra.to(u.arcsec).value,
                                               optimalcube.dec.to(u.arcsec).value, [5, 9])

    dispersionmap_inter, _, _ = Img_interpsmooth(dispersionmap,
                                                 optimalcube.ra.to(u.arcsec).value,
                                                 optimalcube.dec.to(u.arcsec).value, [5, 9])

    optimalimg_inter[snrmap_inter < noise_level] = 0
    velocitymap_inter[snrmap_inter < noise_level] = np.nan
    dispersionmap_inter[snrmap_inter < noise_level] = np.nan

    return optimalimg_inter, velocitymap_inter, \
           dispersionmap_inter, snrmap_inter, x_inter, y_inter

# #定义每一种cube的边界
ly_dataedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':4002*u.AA,'zhi':4040*u.AA}
ly_noiseedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':4050*u.AA,'zhi':4100*u.AA}
ly_continuumedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':4100*u.AA,'zhi':4300*u.AA}
#
lymap=Map(filename='/Users/shiwuzhang/WS/ASTRO/MAMMOTH_KCWI/1441+4003_comb_ss_icubes.fits',
          refpoint=np.array([220.3520886,40.05269183])*u.deg,dataedge=ly_dataedge,
          continuumedge=ly_continuumedge,noiseedge=ly_noiseedge,
          noiselevel=1.8)


lyoptimalcube,lynoisecube,_,lyemissioncube=lymap.optimalmap(smooth=False)
lysnrmap=lymap.snrmap()
lyvelocitymap=lymap.momentmap(order=1,redshift=2.31,restvalue=1215.67*u.AA)
lydispersionmap=lymap.momentmap(order=2,redshift=2.31,restvalue=1215.67*u.AA)


optimalimg,velocitymap,dispersionmap,snrmap,x,y=preprocess(lyoptimalcube,lyvelocitymap,
                                                          lydispersionmap,lysnrmap,3.)
plt.imshow(velocitymap,cmap='Spectral_r')
plt.show()


# heii_dataedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5465*u.AA,'zhi':5495*u.AA}
# heii_noiseedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5400*u.AA,'zhi':5450*u.AA}
# heii_continuumedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5500*u.AA,'zhi':5600*u.AA}
# heiimap=Map(filename='/Users/shiwuzhang/WS/ASTRO/MAMMOTH_KCWI/1441+4003_comb_psfs_icubes.fits',
#             refpoint=np.array([220.3520886, 40.05269183]) * u.deg,dataedge=heii_dataedge,
#             continuumedge=heii_continuumedge,noiseedge=heii_noiseedge,
#             noiselevel=1)
# heiioptimalcube,_,_,heiiemissioncube=heiimap.optimalmap(smooth=False)
# heiisnrmap=heiimap.snrmap()
# heiivelocitymap=heiimap.momentmap(order=1,redshift=2.34,restvalue=1640*u.AA)
# heiidispersionmap=heiimap.momentmap(order=2,redshift=2.34,restvalue=1640*u.AA)
# heiisnrmap[heiisnrmap<3]=0
#
# p_map=np.random.rand(heiisnrmap.shape[0],heiisnrmap.shape[1])
# bool_map=p_map>0.95
# bool_map[20:45,10:]=True
#
# heiisnrmap[~bool_map]=0
#
# optimalimg,velocitymap,dispersionmap,snrmap,x,y=preprocess(heiioptimalcube,heiivelocitymap,
#                                                           heiidispersionmap,heiisnrmap,noise_level=3.)
#
#
# optimalimg[optimalimg==0]=np.power(10,-2.8)
# optimalimg=optimalimg*1e3
# noise=np.random.normal(0,.01,optimalimg.shape)
# optimalimg+=noise


# civ_dataedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5165*u.AA,'zhi':5200*u.AA}
# civ_noiseedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5100*u.AA,'zhi':5150*u.AA}
# civ_continuumedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5000*u.AA,'zhi':5100*u.AA}
#
# civmap=Map(filename='/Users/shiwuzhang/WS/ASTRO/MAMMOTH_KCWI/1441+4003_comb_psfs_icubes.fits',
#            refpoint=np.array([220.3520886, 40.05269183]) * u.deg,dataedge=civ_dataedge,
#            continuumedge=civ_continuumedge,noiseedge=civ_continuumedge,
#            noiselevel=1.5)
#
# civoptimalcube,_,_,civemissioncube=civmap.optimalmap(smooth=True)
# civsnrmap=civmap.snrmap()
# civsnrmap[civsnrmap<3]=0
#
#
# p_map=np.random.rand(civsnrmap.shape[0],civsnrmap.shape[1])
# bool_map=p_map>0.95
# bool_map[20:46,10:-2]=True
#
# civsnrmap[~bool_map]=0
#
# plt.imshow(civsnrmap)
# plt.show()