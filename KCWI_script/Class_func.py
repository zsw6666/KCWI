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
This script define 2 class and 1 functions,
Class:
    1) Cube: inherit from spectral_cube.SpectralCube.
    Implement it with new functions
    2) Map: generate some maps(SNR map, moment map and optimal-extracted map)
    of datacube based on Cube
    3) Img_interpsmooth: interpolate and smooth the images.
'''




import numpy as np
from astropy import units as u
from astropy.convolution import convolve,Gaussian2DKernel
from spectral_cube import SpectralCube
from scipy import ndimage



class Cube(SpectralCube):

    def __init__(self,refpoint=[0*u.deg,0*u.deg],slices=None,
                 optimal_img=None,continuum=None,**krg):
        super().__init__(**krg)
        self.refpoint=refpoint
        self.dec = np.mean(self.spatial_coordinate_map[0] - self.refpoint[1], axis=0)
        self.ra = np.mean(self.spatial_coordinate_map[1] - self.refpoint[0], axis=1)
        self.continuum=continuum
        self.slices=slices
        self.optimal_img=optimal_img

    def emissionline_subcube(self,mean_std):
        '''
        extract the best subcube which fully include the
        emission line but doesn't have continuum component.
        Idea is simple, calculate standard deviation of each
        slice, if there's emission line, standard deviation would
        be large than other slices.
        '''

        #select slices with emission line component
        slices_std=self.noise()
        slices_bool=(slices_std>mean_std)
        presubcube_axis = self.spectral_axis[slices_bool]

        #choose the longest consecutive array and define it
        #as the emission-line interval
        step_size=(self.spectral_axis[1]-self.spectral_axis[0])
        consecutive_arr=np.split(presubcube_axis,np.where(np.diff(presubcube_axis)>1.5*step_size)[0]+1)
        consecutive_arr_size=np.array([len(i) for i in consecutive_arr])
        longest_arr=consecutive_arr[np.where(consecutive_arr_size==consecutive_arr_size.max())[0][0]]

        emissionline_cube=self.subcube(zlo=longest_arr[0],zhi=longest_arr[-1])
        return emissionline_cube

    def estimate_continuum(self):
        '''
        if this cube is a continuum cube,
        calculate the continuum image with
        this cube
        '''
        self.continuum=self.mean(axis=0)

        return None

    def mask_generate(self,noise,step=5):
        '''
        generate spatial and spectral mask to
        extract the signal beyond noise*noise_level
        from each slice or each pixel. noise is an
        1D or 2D array, it should have same shape with
        spectral axis or spatial axis of self.cube
        for each slice, there's a noise
        '''

        maskcube=np.zeros_like(self.hdu.data)
        slicecube=np.add.reduceat(self.hdu.data,range(0,self.hdu.data.shape[0],step),axis=0)
        weight=np.add.reduceat(np.ones(self.hdu.data.shape[0]),range(0,self.hdu.data.shape[0],step))
        noise=noise*weight

        noise=noise.reshape(slicecube.shape[0],1,1)
        slicecube_bool=(slicecube-noise>=0)
        i,j=0,0
        while i<maskcube.shape[0]:
            if weight[j]>0:
                maskcube[i][slicecube_bool[j]]=1
                weight[j]-=1
                i+=1
            else:
                j+=1

        return maskcube

    def max_emission_extraction(self,maskcube):
        '''
        extract the emission region as large as possible.
        The idea is to select
        '''

        optimal_datacube=self.hdu.data*maskcube
        mask=np.sum(maskcube,axis=0)
        mask[mask==0]=1.
        self.optimal_img=np.sum(optimal_datacube,axis=0)/mask
        optimalcube=Cube(data=optimal_datacube,wcs=self.wcs,
                         refpoint=self.refpoint)
        optimalcube.optimal_img=self.optimal_img
        return optimalcube

    def noise(self):
        '''
        if cube is a noise cube, calculate
        standard deviation for each slice and
        return a noise array, this array should be
        a 1D array with the same length with noise cube
        '''

        noise=np.zeros(self.hdu.data.shape[0])
        for i in range(self.hdu.data.shape[0]):
            std=np.std(self.hdu.data[i])
            noise[i]=std
        return noise

    @classmethod
    def read(cls,filename,**krg):
        '''
        class function, read the datacube from .fits
        '''
        cube=super().read(filename=filename)
        cube=Cube(data=cube.hdu.data,wcs=cube.wcs,**krg)
        return cube

    def slice_generate(self,step):
        '''
        choose a interval and step, return a series of
        slices in the interval with given step
        '''

        '''
        extract slices in the loop
        '''
        IMG,SLICE_loc=[],[]
        i=0
        while i<self.shape[0]-step:
            img=np.sum(self[i:i+step,:,:].hdu.data,axis=0)
            IMG.append(img)
            location=np.mean(np.mean(self.spectral_axis[i:i+step]).value)
            SLICE_loc.append(location)
            i += step

        #location of the slices
        self.slices=(IMG,SLICE_loc)
        return None

    def substract_continuum(self,continuum_img):
        '''
        subtract continuum component from self.
        continuum_img should have the same spatial
        shape like self
        '''
        continuum_img=continuum_img.reshape(1,continuum_img.shape[0],
                                            continuum_img.shape[1])
        continuum_sub=self-continuum_img

        return continuum_sub

    def snrmap(self,noise,maskcube):
        '''
        Calculate the SNR map for
        optimal-extracted image
        '''
        if self.optimal_img is not None:
            mask=np.sum(maskcube,axis=0)
            snrmap=(self.optimal_img-0.7*noise)*np.sqrt(mask)/noise
        else:
            print("There's no optimal-extracted image for cube!")
        return snrmap

    def subcube(self, xlo='min', xhi='max', ylo='min', yhi='max', zlo='min',
                zhi='max', rest_value=None, **krg):

        presub_cube=super().subcube(xlo,xhi,ylo,yhi,zlo,zhi,rest_value)
        sub_cube=Cube(data=presub_cube.hdu.data,wcs=presub_cube.wcs,refpoint=self.refpoint)

        return sub_cube

class Map:

    def __init__(self,filename,refpoint,
                 dataedge,continuumedge,noiseedge,noiselevel):

        self.cube=Cube.read(filename,refpoint=refpoint)
        self.refpoint=refpoint
        self.dataedge=dataedge
        self.continuumedge = continuumedge
        self.noiseedge = noiseedge
        self.noiselevel=noiselevel

    def _basecube_generate(self):
        '''
        This function is the inner function
        which generate three base cubes(datacube,
        noisecube and continuumcube) for map
        '''

        datacube=self.cube.subcube(refpoint=self.refpoint,**self.dataedge)
        noisecube=self.cube.subcube(refpoint=self.refpoint,**self.noiseedge)
        continuumcube = self.cube.subcube(refpoint=self.refpoint,**self.continuumedge)


        return datacube,noisecube,continuumcube

    def optimalmap(self):
        '''
        generate the best extracted psudo-nb image
        '''
        datacube,noisecube,continuumcube=self._basecube_generate()
        emissioncube = datacube.emissionline_subcube(noisecube.noise().mean()*self.noiselevel['emission'])
        continuumcube.estimate_continuum()
        maskcube=emissioncube.mask_generate(self.noiselevel['mask']*noisecube.noise().mean())
        optimalcube=emissioncube.max_emission_extraction(maskcube)
        optimalcube.optimal_img=optimalcube.optimal_img+noisecube.noise().mean()*0.7-continuumcube.continuum
        return optimalcube,noisecube,maskcube,emissioncube

    def snrmap(self):
        '''
        calculate the signal-to-noise map for
        optimal-extracted image
        '''
        optimalcube,noisecube,maskcube,_=self.optimalmap()
        snrmap = optimalcube.snrmap(noisecube.noise().mean(), maskcube)
        return snrmap

    def momentmap(self,order,redshift,restvalue):
        '''
        calculate the first or second order of moment
        map
        '''
        optimalcube, noisecube, maskcube,_ = self.optimalmap()
        optimalcube=optimalcube.with_spectral_unit(u.km/u.s,
                                                   velocity_convention='relativistic',
                                                   rest_value=(1+redshift)*restvalue)
        momentmap=optimalcube.moment(order=order,axis=0)

        momentmap=np.power(momentmap,1/order)

        return momentmap.to(u.km/u.s).value

def Img_interpsmooth(img,x,y,n_inter):
    '''
    interpolate the image with given interpolate
    number n_inter, this should be an 1D array with
    two cell corresponding to x and y respectively.
    Meanwhile, interpolate x, y
    '''

    kernel=Gaussian2DKernel(x_stddev=1,y_stddev=1,x_size=1,y_size=3)
    img_smooth=convolve(img,kernel)
    img_smooth_inter=ndimage.zoom(img_smooth,n_inter)
    x_inter=ndimage.zoom(x,n_inter[0])
    y_inter=ndimage.zoom(y,n_inter[1])

    return img_smooth_inter,y_inter,x_inter
