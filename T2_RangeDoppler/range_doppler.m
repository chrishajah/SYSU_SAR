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
tpd = 30e-6;                       % Pulse width 30us

% Verify that the range resolution is as expected
% 验证距离向分辨率
rngRes1 = bw2rangeres(bw_10x10);
rngRes2 = bw2rangeres(bw_1x1);

% 定义方位向分辨率
aziRes1 = rngRes1;
aziRes2 = rngRes2;

% 根据方位向分辨率计算得到合成孔径
% Antenna properties
La1 = 2*aziRes1;
La2 = 2*aziRes2;


%%
%Radar Platform Parameters
%场景参数：场景中心距离雷达为10Km，方位×距离——120m×520m，雷达高度为0。

H = 0;                          %雷达平台高度  0
theta_a = 0;                    %雷达方位斜视角 0度
Vr = 200;                       %雷达平台的运动速度200m/s
R_eta_c = 10e3;                 %场景中心距离距离10km


%目标参数：25个RCS等于1的理想点目标，均匀分布在以场景中心为中心的方位-50m～+50m、距离-250m～+250m的范围内。
%-------------方位向120m---------------]
%   -50   -25    0     25    50
%
%
%
%    x     x     x     x     x   50
%
%
%    x     x     x     x     x   25      距
% 
%
%    x     x     x     x     x   0       离 (520m)
%
%                                   
%    x     x     x     x     x   -25     向
%
%
%    x     x     x     x     x   -50
%
%
%
%
%-------------方位向120m--------------


targets_pos = [];

