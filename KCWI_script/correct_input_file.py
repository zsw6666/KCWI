# !/usr/bin/python3
# -*- coding: utf-8 -*-
# @Time      :   2020/7/6 1:16 PM
# @Author    :   Shiwu
# @Site      :   
# @File      :   correct_input_file.py
# @Desc      :   
# @Software  :   PyCharm
# @license   :   Copyright(C), C
# @Contact   :   shiwuzhang95@gmail.com


import os
import numpy as np

dic='/Users/shiwuzhang/WS/ASTRO/MAMMOTH_KCWI/models/model_v2/photoion/'
file_list=[file for file in os.listdir(dic) if 'nlr' in file]

for file in file_list:

    with open(dic+file,'r+') as f:
        new_f=f.readlines()
        f.seek(0)
        for line in new_f:
            if 'luminosity' not in line:
                if 'agn' in line:
                    f.write(line[:9]+'\n')
                else:
                    f.write(line)
        f.truncate()

