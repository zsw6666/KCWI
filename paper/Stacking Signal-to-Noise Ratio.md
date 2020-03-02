# Stacking Signal-to-Noise Ratio

It really confused me before why stacking can improve the signal-to-noise-ratio, here I will talk about this problem.

If there's only one slice of image, the routine to calculate SNR for each pixel is given:

$SNR=\frac{S_{ij}}{\sigma}$

$S_{ij}$ is the flux captured by pixel (i,j), $\sigma$ is the background standard deviation which is given by $\sigma=\frac{\sum_{i,j}n_{ij}}{N_{b}}$ , $N_{b}$ is the number of pixels with no signal and  $n_{ij}$ is the value of pixel (i,j).

This is the method to calculate the SNR for a single piece of image. However, the situation of observation is that one pointing of target is always divided into serveral parts which can make the observation flexible. Then All individual image would be stacked together to produce the final image. This process can improve the SNR of our sources.  The Following calculation gives the explaination.

<img src="/Users/shiwuzhang/WS/ASTRO/MAMMOTH_KCWI/draft_code/paper/SNR.jpg" alt="img" style="zoom:30%;" />





