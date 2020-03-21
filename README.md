# Synthetic-Aperture-Radar-SAR-Imaging-algorithms——————chirp scaling algorithm

1）Final_cs_broadside.m : 
working conditions：height of flying platform is 0；broadside mode；slant plane imaging
（不考虑高度、正侧视、斜平面成像、cs算法）

2）Final_cs_squint5.m : 
working conditions：height of flying platform is 0；squint angle is 5 degree；slant plane imaging.
The imaging results are poor for the range sampling rate (Fs = 200e6) is too low, which leads to the folding of two-dimension spectrum.
不考虑高度、斜视 5 度、斜平面cs成像算法，效果不好；因为二维频谱折叠，需提高距离向采样率 Fs

Improved versions for squint mode imaging:
3)Final_cs_squint5_new.m : 
working conditions：height of flying platform is 0；squint angle is 5 degree；slant plane imaging; Fs = 400e6.
不考虑高度、斜视 5 度、斜平面cs成像算法，效果好

4)Final_cs_squint5_new_1.m : 
working conditions：height of flying platform is 0；squint angle is 5 degree；slant plane imaging; Fs = 200e6.
Range upsampling is performed to avoid two-spectrum folding.
不考虑高度、斜视 5 度、斜平面cs成像算法，效果好
