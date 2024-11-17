%SAR Range-Doppler Algorithm Simulation
%ChrisL 2024 @ SYSU
clear;clc;

%%
%Define SAR System Parameters
%LFM信号参数：中心频率10.0GHz，脉冲宽度30us。
%SAR的方位斜视角等于0，即正侧视。
%天线加权形式：矩形加权或sinc加权。
%分辨率有两种情况：（1）10m×10m；（2）1m×1m。

freq = 10e9;                       % Carrier frequency 10GHz
[lambda,c] = freq2wavelen(freq);   % Wavelength (m) 
bw_10x10 = 14.985e6;               % Signal bandwidth for 10x10 range resolution -- 14.985Mhz
bw_1x1 = 14.985e7;                 % Signal bandwidth for 1x1 range resolution -- 149.85Mhz
fs = 60e6;                         % Sampling frequency (Hz)
tpd = 30e-6;                       % Pulse width 30us
Kr_10x10 = bw_10x10 * tpd;
Kr_1x1 = bw_1x1 * tpd;
% Verify that the range resolution is as expected
rngRes1 = bw2rangeres(bw_10x10);
rngRes2 = bw2rangeres(bw_1x1);

alpha_s_r = 1.2;             % 距离向过采样因子
alpha_s_a = 1.2;             % 方位向过采样因子

% Antenna properties
apertureLength_10x10 = 30;         % Aperture length for 10x10 resolution(m) 
apertureLength_1x1 = 300;          % Aperture length for 1x1 resolution(m) 
sqa = 0;                           % Squint angle (deg) SAR的方位斜视角等于0，即正侧视



%%