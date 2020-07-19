# !/usr/bin/python3
# -*- coding: utf-8 -*-
# @Time      :   2020/6/9 8:55 AM
# @Author    :   Shiwu
# @Site      :   
# @File      :   IO_mysql.py
# @Desc      :   
# @Software  :   PyCharm
# @license   :   Copyright(C), C
# @Contact   :   shiwuzhang95@gmail.com

import pandas as pd
import pymysql
import numpy as np

#information for login to the mysql.
host='3mdb.astro.unam.mx'
user='OVN_user'
passwd='oiii5007'
port='3306'
DB_s='3MdBs'
DB_p_17='3MdB_17'

#command to extract the data
sql_command_origin='''
select shck_vel as v, Z as z, 
CIV_1548 as CIV, HeII_1640 as HeII, HI_1216 as ly, HI_6563 as h_alpha,
preshck_dens as density, preshck_temp as T, distance as D,
cut_off_temp as T_c,mag_fld as B,time as t,shock_params.ref,
FHI,FHII,FHeI,FHeII,FHeIII,ProjectID,script
from shock_params
inner join emis_UVC on emis_UVC.ModelID=shock_params.ModelID 
inner join emis_UVB on emis_UVB.ModelID=shock_params.ModelID 
inner join emis_VI on emis_VI.ModelID=shock_params.ModelID 
inner join abundances on abundances.AbundID=shock_params.AbundID 
where emis_UVC.model_type='shock'
and emis_UVB.model_type='shock' 
and emis_VI.model_type='shock'
and Z=
'''

#the metallicity
savdic='/work/zsw666/MAMMOTH_KCWI/3MdB_mysql/'
metal_Z=pd.read_csv(savdic+'metal_z.csv')

#extract data for all metallicity
savdic='/work/zsw666/MAMMOTH_KCWI/3MdB_mysql/model/SP'
print('connecting.....')
db = pymysql.connect(host=host, db=DB_s, user=user, passwd=passwd,
                     connect_timeout=2400, max_allowed_packet=16777216*20)
k=1
for i in metal_Z['Z']:
    sql_command=sql_command_origin[:-1]+'%f'%i
    print('export data.....%d'%k)
    res = pd.read_sql(sql_command, con=db)
    file_name='SP_Z_%.3f.csv'%np.around(i/0.02,decimals=3)
    res.to_csv(savdic+file_name,encoding='utf-8', index=False)
    k+=1
db.close()
