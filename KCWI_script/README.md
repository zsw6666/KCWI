## ASTRO_script

**These script do data analysis for KCWI data cube.**

here're brief introduction to these script.

1. Class_func:
   * Includes 2 calss and 1 function.
     * Class Cube: inherit from spectral_cube.SpectralCube. Implement it with some new functions
     * Class Map: produce maps including SNR map, moment map and optimal-extracted map based on Cube.
     * Function Img_interpsmooth: interpolate and smmoth images.
2. Map:
   * generate SNR map,  moment map, optimal-extracted map with KCWI data cube.
3. Contour_generator:
   * produce images for futhermore contour overlay on HST images.
4. Contourplot:
   * plot HST image and overlay contour on it.
5. Emissionflux:
   * calculate the luminosity for ly$\alpha$ HeII and CIV
6. Emissionfit:
   * Fit the spectral extracted from different part of ly$\alpha$ nebulae.
7. Map_generator:
   * Produce interpolated and smoothed map from Map.py
8. Mapplot:
   * Plot the interpolated and smoothed map from Map_generator.
9. Slicesplot:
   * Plot slices of ly$\alpha$ emission.
10. Specmap:
    * extract spectra from different part of ly$\alpha$ nebulae.
11. Specmapplot:
    * plot spectra from spec map.
12. Spectralplot:
    * plot spectra extracted center on source-B.

