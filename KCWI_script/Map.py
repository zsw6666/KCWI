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
from Class_func import Map


#定义每一种cube的边界
ly_dataedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':4002*u.AA,'zhi':4040*u.AA}
ly_noiseedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':4050*u.AA,'zhi':4100*u.AA}
ly_continuumedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':4100*u.AA,'zhi':4300*u.AA}
heii_dataedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5465*u.AA,'zhi':5495*u.AA}
heii_noiseedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5400*u.AA,'zhi':5450*u.AA}
heii_continuumedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5500*u.AA,'zhi':5600*u.AA}
civ_dataedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5165*u.AA,'zhi':5200*u.AA}
civ_noiseedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5100*u.AA,'zhi':5150*u.AA}
civ_continuumedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':5000*u.AA,'zhi':5100*u.AA}

lymap=Map(filename='/Users/shiwuzhang/WS/ASTRO/MAMMOTH_KCWI/1441+4003_comb_ss_icubes.fits',
          refpoint=np.array([220.3520886,40.05269183])*u.deg,dataedge=ly_dataedge,
          continuumedge=ly_continuumedge,noiseedge=ly_noiseedge,
          noiselevel=.2)

heiimap=Map(filename='/Users/shiwuzhang/WS/ASTRO/MAMMOTH_KCWI/1441+4003_comb_psfs_icubes.fits',
            refpoint=np.array([220.3520886, 40.05269183]) * u.deg,dataedge=heii_dataedge,
            continuumedge=heii_continuumedge,noiseedge=heii_noiseedge,
            noiselevel=.2)

civmap=Map(filename='/Users/shiwuzhang/WS/ASTRO/MAMMOTH_KCWI/1441+4003_comb_psfs_icubes.fits',
           refpoint=np.array([220.3520886, 40.05269183]) * u.deg,dataedge=civ_dataedge,
           continuumedge=civ_continuumedge,noiseedge=civ_continuumedge,
           noiselevel=.2)


lyoptimalcube,lynoisecube,_,lyemissioncube=lymap.optimalmap()
lysnrmap=lymap.snrmap()
lyvelocitymap=lymap.momentmap(order=1,redshift=2.31,restvalue=1215.67*u.AA)
lydispersionmap=lymap.momentmap(order=2,redshift=2.31,restvalue=1215.67*u.AA)


heiioptimalcube,_,_,heiiemissioncube=heiimap.optimalmap()
heiisnrmap=heiimap.snrmap()
heiivelocitymap=heiimap.momentmap(order=1,redshift=2.34,restvalue=1640*u.AA)
heiidispersionmap=heiimap.momentmap(order=2,redshift=2.34,restvalue=1640*u.AA)


civoptimalcube,_,_,civemissioncube=civmap.optimalmap()
civsnrmap=civmap.snrmap()
civvelocitymap=civmap.momentmap(order=1,redshift=2.34385,restvalue=1549*u.AA)
civdispersionmap=civmap.momentmap(order=2,redshift=2.34385,restvalue=1549*u.AA)
