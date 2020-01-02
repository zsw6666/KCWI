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


import numpy as np
from astropy import units as u
from Class_func import Map
# s

ly_dataedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':4000*u.AA,'zhi':4100*u.AA}
ly_noiseedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':4050*u.AA,'zhi':4100*u.AA}
ly_continuumedge={'xlo':2,'xhi':22,'ylo':4,'yhi':65,'zlo':4100*u.AA,'zhi':4300*u.AA}
heii_dataedge={'xlo':10,'xhi':21,'ylo':28,'yhi':59,'zlo':5450*u.AA,'zhi':5500*u.AA}
heii_noiseedge={'xlo':10,'xhi':21,'ylo':28,'yhi':59,'zlo':5400*u.AA,'zhi':5450*u.AA}
heii_continuumedge={'xlo':10,'xhi':21,'ylo':28,'yhi':59,'zlo':5500*u.AA,'zhi':5600*u.AA}
civ_dataedge={'xlo':10,'xhi':21,'ylo':30,'yhi':59,'zlo':5100*u.AA,'zhi':5215*u.AA}
civ_noiseedge={'xlo':10,'xhi':21,'ylo':30,'yhi':59,'zlo':5100*u.AA,'zhi':5150*u.AA}
civ_continuumedge={'xlo':10,'xhi':21,'ylo':30,'yhi':59,'zlo':5000*u.AA,'zhi':5100*u.AA}

lymap=Map(filename='/Users/shiwuzhang/WS/ASTRO/MAMMOTH_KCWI/1441+4003_comb_ss_icubes.fits',
          refpoint=np.array([220.3520886,40.05269183])*u.deg,dataedge=ly_dataedge,
          continuumedge=ly_continuumedge,noiseedge=ly_noiseedge,
          noiselevel={'emission':1,'mask':1.})

heiimap=Map(filename='/Users/shiwuzhang/WS/ASTRO/MAMMOTH_KCWI/1441+4003_comb_psfs_icubes.fits',
            refpoint=np.array([220.3520886, 40.05269183]) * u.deg,dataedge=heii_dataedge,
            continuumedge=heii_continuumedge,noiseedge=heii_noiseedge,
            noiselevel={'emission':0.98,'mask':1.})

civmap=Map(filename='/Users/shiwuzhang/WS/ASTRO/MAMMOTH_KCWI/1441+4003_comb_psfs_icubes.fits',
           refpoint=np.array([220.3520886, 40.05269183]) * u.deg,dataedge=civ_dataedge,
           continuumedge=civ_continuumedge,noiseedge=civ_continuumedge,
           noiselevel={'emission': 1, 'mask': 1.})


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
