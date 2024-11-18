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
Kr = bw_1x1 * tpd;

% Verify that the range resolution is as expected
% 验证距离向分辨率
rngRes1 = bw2rangeres(bw_10x10);
rngRes2 = bw2rangeres(bw_1x1);

% 定义方位向分辨率
aziRes1 = rngRes1;
aziRes2 = rngRes2;

% 根据方位向分辨率计算得到真实孔径
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

%根据真实孔径计算合成孔径大小
Ls = 0.886*lambda*R_eta_c / La2;
%计算多普勒频宽
bw_dop = 0.886*2*Vr*cos(theta_a)/La2;
f_eta_c = 2 * Vr * sin(theta_a) / lambda; % 多普勒中心频率


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


% 过采样率设置为1.2
alpha_o = 1.2;

% 方位向慢时间采样率
f_azi = alpha_o*bw_dop; 

% 距离向快时间采样率
f_rng = alpha_o*bw_1x1;


Nrg = 8192;
Naz = 128;
%定义快慢时间时域变量
eta = -60/Vr:1/f_azi:60/Vr;
tau = -15e-6:1/f_rng:15e-6;


%定义频域变量
f_tau = (-f_rng / 2: f_rng/Nrg :f_rng / 2 - f_rng/Nrg); % 距离频率变量
f_tau=f_tau-(round(f_tau/f_rng))/f_rng;%混叠方程
f_eta = f_eta_c + (-f_azi / 2:f_azi/Naz:f_azi / 2 - f_azi/Naz); % 方位频率变量
f_eta=f_eta-(round(f_eta/f_azi))/f_azi;

%%
%信号建模(先对单个点目标0,0回波信号建模)

[Tau,Eta] = meshgrid(tau,eta);
[F_tau, F_eta] = meshgrid(f_tau, f_eta); 


R_Eta = sqrt(R_eta_c^2 + Vr^2 * Eta.^2);

Wr = (abs(Tau-2*R_Eta/c) <= tpd / 2);
Wa = sinc((La2*atan(Vr*(Eta - 0)./R_eta_c)/lambda).^2);

S_0 = Wa.*exp(-1j*4*pi*f_rng .*R_Eta./c).*exp(1j*pi*Kr.*(Tau - 2.*R_Eta./c).^2);

figure();
imagesc(abs(S_0));

%距离压缩滤波器
Hrg = (abs(F_tau) <= bw_1x1 / 2) .* exp(+1j*pi*F_tau.^2/Kr);%滤波器
S_f = fft(S_0,8192,2);

S_rngcomp = fftshift(ifft(Hrg .* S_f));
S_f_rngcomp = Hrg .* S_f;


figure();
imagesc(abs(S_rngcomp));
