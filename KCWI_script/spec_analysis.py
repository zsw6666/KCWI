# !/usr/bin/python3
# -*- coding: utf-8 -*-
# @Time      :   2020/7/9 1:28 PM
# @Author    :   Shiwu
# @Site      :   
# @File      :   spec_analysis.py
# @Desc      :   
# @Software  :   PyCharm
# @license   :   Copyright(C), C
# @Contact   :   shiwuzhang95@gmail.com


import numpy as np
from scipy.signal import find_peaks
from lmfit.models import GaussianModel,ConstantModel

from spec_extraction import spec_ly_lib,x_window_ly,y_window_ly,r_ly \
, spec_heii_lib,x_window_heii,y_window_heii,r_heii \
, spec_civ_lib,x_window_civ,y_window_civ,r_civ \
, spec_ly_lib2,x_window_ly2,y_window_ly2,r_ly2 \
, v_ly,v_heii,v_civ,snr_lib,v_lib,ly_noise,heii_noise,civ_noise

spec_ly_lib=spec_ly_lib*1e3
spec_heii_lib=spec_heii_lib*1e3
spec_civ_lib=spec_civ_lib*1e3

ly_noise=ly_noise*1e3/5.
heii_noise=heii_noise*1e3/5.
civ_noise=civ_noise*1e3/5.

class kinematics:
    def __init__(self, speclib, x):
        self.y = speclib
        self.x = x

    def velocity(self, y_index):
        y = self.y[y_index] / self.y[y_index].max()
        v_aver = (self.x * y).sum() / y.sum()
        return v_aver

    def dispersion(self, y_index):
        v_aver = self.velocity(y_index)
        y = self.y[y_index] / self.y[y_index].max()
        variance = (((self.x - v_aver) ** 2) * y).sum() / y.sum()
        sigma = np.sqrt(np.abs(variance))

        return sigma

    def window(self):
        v_lib, sigma_lib = [], []
        for i in range(self.y.shape[0]):
            v_lib.append(self.velocity(i))
            sigma_lib.append(self.dispersion(i))

        return np.array(v_lib), np.array(sigma_lib)

class fittor:

    def __init__(self, speclib, lamda):
        self.speclib = speclib
        self.lamda = lamda

    def gaussian_fit(self, x, y, paras=None, method=None):
        g1model = GaussianModel(prefix='g1_')
        g2model = GaussianModel(prefix='g2_')
        g3model = GaussianModel(prefix='g3_')
        cmodel = ConstantModel()

        paras_fit = g1model.guess(data=y, x=x)
        paras_fit.update(g2model.make_params())
        paras_fit.update(g3model.make_params())

        paras_fit['g1_amplitude'].set(min=0.)
        paras_fit['g2_amplitude'].set(min=0.)
        paras_fit['g3_amplitude'].set(min=0.)

        paras_fit['g1_center'].set(min=paras['g1_center'] - 300.,
                                   value=paras['g1_center'],
                                   max=paras['g1_center'] + 300.)
        paras_fit['g2_center'].set(min=paras['g2_center'] - 300.,
                                   value=paras['g2_center'],
                                   max=paras['g2_center'] + 300.)
        paras_fit['g3_center'].set(min=paras['g3_center'] - 300.,
                                   value=paras['g3_center'],
                                   max=paras['g3_center'] + 300.)

        paras_fit['g1_sigma'].set(min=100, value=paras['g1_sigma'],
                                  max=paras['g1_sigma'] + 50.)
        paras_fit['g2_sigma'].set(min=100, value=paras['g2_sigma'],
                                  max=paras['g2_sigma'] + 50.)
        paras_fit['g3_sigma'].set(min=100, value=paras['g3_sigma'],
                                  max=paras['g3_sigma'] + 50.)

        paras_fit.update(cmodel.make_params())

        model = g1model + g2model + g3model + cmodel
        result = model.fit(y, x=x, params=paras_fit, mothod=method)

        yfit = result.best_fit
        y_para = result.best_values

        residual = y - yfit

        return yfit, y_para, result, residual

    def _init_parameter(self, x, y):

        peaks, properties = find_peaks(y, prominence=1., width=3, distance=3)
        v = np.zeros(3)
        peaks = list(peaks)
        for i in range(3):
            if len(peaks) > 0:
                v[i] = x[peaks.pop()]
            else:
                break
        init_para = {'g1_center': v[0], 'g2_center': v[1], 'g3_center': v[2],
                     'g1_sigma': 200., 'g2_sigma': 200., 'g3_sigma': 200.}

        return init_para

    def best_fit(self, x, y, criteria, method, inter=50):

        sigma_per, k = 0., 0
        paras = self._init_parameter(x, y)
        while sigma_per < 0.6:
            yfit, ypara, out, residual = self.gaussian_fit(x, y, paras=paras, method=method)
            sigma_per = residual[(residual >= -criteria) & (residual <= criteria)].size / residual.size
            paras = ypara
            k += 1
            if k >= inter:
                break

        return yfit, ypara, out, residual, sigma_per

    def lib_fit(self, method, criteria=0.55):
        fitlib = []
        paralib = []
        outlib = []
        residuallib = []
        sigma_per_lib = []
        for i in range(self.speclib.shape[0]):
            yfit, para, out, residual, sigma_per = self.best_fit(self.lamda, self.speclib[i, :], criteria, method)
            fitlib.append(yfit)
            paralib.append(para)
            outlib.append(out)
            residuallib.append(residual)
            sigma_per_lib.append(sigma_per)

        return np.array(fitlib), paralib, np.array(residuallib), outlib, np.array(sigma_per_lib)

ly_kinematic=kinematics(spec_ly_lib,v_ly)
v_window_ly,sig_window_ly=ly_kinematic.window()

heii_kinematic=kinematics(spec_heii_lib,v_heii)
v_window_heii,sig_window_heii=heii_kinematic.window()

civ_kinematic=kinematics(spec_civ_lib,v_civ)
v_window_civ,sig_window_civ=civ_kinematic.window()


ly_fittor=fittor(spec_ly_lib,v_ly)
fitted_ly,para_ly,residual_ly,out_ly,sigma_per=ly_fittor.lib_fit(method='lbfgsb')

heii_fittor=fittor(spec_heii_lib,v_heii)
fitted_heii,para_heii,residual_heii,out_heii,sigma_per=heii_fittor.lib_fit(method='lbfgsb',criteria=0.25)

civ_fittor=fittor(spec_civ_lib,v_civ)
fitted_civ,para_civ,residual_civ,out_civ,sigma_per=civ_fittor.lib_fit(method='lbfgsb',criteria=.28)
