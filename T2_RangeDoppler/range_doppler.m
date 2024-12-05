%SAR Range-Doppler Algorithm Simulation
%ChrisL 2024 @ SYSU
clear;clc;close all;

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
% 雷达3dB波束
beta_bw = Ls/R_eta_c;
%计算多普勒频宽
bw_dop = 0.886*2*Vr*cos(theta_a)/La2;
f_eta_c = 2 * Vr * sin(theta_a) / lambda; % 多普勒中心频率

% 过采样率设置为1.2
alpha_o = 1.2;

% 方位向慢时间采样率
f_azi = alpha_o*bw_dop; 

% 距离向快时间采样率
f_rng = alpha_o*bw_1x1;

t_eta_c = -R_eta_c * sin(theta_a) / Vr; % 波束中心穿越时刻

% 定义快慢时间时域变量

% 定义方位和距离向频域FFT点数
Nrg = 8192;
Naz = 256;

Trg = Nrg / f_rng;Taz = Naz / f_azi;
%距离向/方位向采样时间间隔
Gap_t_tau = 1 / f_rng;Gap_t_eta = 1 / f_azi;
%距离向/方位向采样频率间隔
Gap_f_tau = f_rng / Nrg;Gap_f_eta = f_azi / Naz;


% eta为慢时间变量，对应方位向 长度：256
t_eta = t_eta_c + (-Taz / 2:Gap_t_eta:Taz / 2 - Gap_t_eta); % 方位时间变量
% tau为快时间变量，对应距离向 长度：8192
t_tau = 2 * R_eta_c / c + (-Trg / 2:Gap_t_tau:Trg / 2 - Gap_t_tau); % 距离时间变量


% 定义频域变量
% 距离频率变量 长度：8192
f_tau = (-f_rng / 2: f_rng/Nrg :f_rng / 2 - f_rng/Nrg);       
f_tau=f_tau-(round(f_tau/f_rng))/f_rng;                         % 混叠方程

% 方位频率变量 长度：256
f_eta = f_eta_c + (-f_azi / 2:f_azi/Naz:f_azi / 2 - f_azi/Naz); 
f_eta=f_eta-(round(f_eta/f_azi))/f_azi;                         % 混叠方程


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

% 目标个数25
target_Num = 25;

% 目标方位向坐标
target_Pos_azimuth = [-50,-25,0,25,50,-50,-25,0,25,50,-50,-25,0,25,50,-50,-25,0,25,50,-50,-25,0,25,50];
target_Pos_range = [50,25,0,-25,-50,50,25,0,-25,-50,50,25,0,-25,-50,50,25,0,-25,-50,50,25,0,-25,-50];
target_R0 =  sqrt(target_Pos_range .^2+R_eta_c^2); %对每个目标计算瞬时斜距
Ta = 0.886 * lambda * R_eta_c / (La2 * Vr * cos(theta_a)); %目标照射时间


%%
% 信号建模

% 设置距离时域-方位时域二维网络坐标
[T_tau,T_eta] = meshgrid(t_tau,t_eta);

% 设置频率时域-方位频域二维网络坐标
[F_tau, F_eta] = meshgrid(f_tau, f_eta); 

% 定义回波时域函数矩阵
S_echo = zeros(Naz, Nrg);

for ii=1:target_Num
    %eta_c_target(ii) = (target_Pos_azimuth(ii)-target_Pos_range(ii)*tan(theta_a))/Vr;% 目标i的波束中心穿越时刻
    eta_c_target(ii) = (target_Pos_azimuth(ii)-target_R0(ii) * tan(theta_a)) / Vr; %波束中心穿越时刻
    R_eta(:,:,ii) = sqrt((target_R0(ii))^2+(Vr^2)*((T_eta - target_Pos_azimuth(ii) / Vr).^2));
    % R_eta(:,:,ii) = sqrt(R_eta_c^2+(target_Pos_range(ii)-(-R_eta_c*sin(theta_a)))^2+(target_Pos_azimuth(ii)-T_eta*Vr).^2);
    %R_eta(:,:,ii) =  target_R0(ii) + Vr^2 .* T_eta .^2 /2*R_eta_c;
    % 距离向包络，即距离窗
    Wr = abs(T_tau - 2 * R_eta(:,:,ii)/c) <=  tpd/2 ;
    % 方位向包络，采用sinc包络
    %Wa = (La2*atan(Vr*(T_eta - eta_c_target(ii))./(R_eta_c * sin(theta_a) + target_Pos_range(ii))/lambda).^2)<=Ta/2;
    Wa = abs(T_eta-eta_c_target(ii)) <= Ta / 2;
    %theta = atan( Vr.*(t_tau-eta_c_target(ii).*ones(Naz,Nrg))/target_Pos_range(ii) );
    %Wa = (sinc(0.886.*theta./beta_bw)).^2; 
    % 回波相位项
    Phase = exp(-1j*4*pi*freq*R_eta(:,:,ii)/c) .* exp(+1j*pi*Kr*(T_tau - 2 * R_eta(:,:,ii) / c).^2);
    S_echo_Target = Wr  .* Wa .* Phase;
    S_echo = S_echo + S_echo_Target;
    %%%%%%%%%%%%%%格式化输出%%%%%%%%%%%%%%%%%%
    fprintf("当前目标:%d,坐标:(%.2f,%.2f),\n最近斜距R0=%.2f,波束中心穿越时刻=%.4f\n", ii, target_Pos_range(ii), target_Pos_azimuth(ii),target_R0(ii),eta_c_target(ii));
    %Time_tau_Point =  round(((2*target_R0)/(c*cos(theta_a))-t_tau(1))/Gap_t_tau);%目标的距离向标准坐标(距离门)
    %Time_eta_Point = Naz / 2 + (target_Pos_azimuth(ii) / Vr) /Gap_t_eta;
    %fprintf("徙动校正前,点数坐标应为%d行(方),%d列(距)\n\n",round(Time_eta_Point(ii)),Time_tau_Point);
end

S_echo = S_echo .* exp(-2j*pi*f_eta_c*T_eta); %校正时间轴


%% 距离压缩

Hf = (abs(F_tau) <= bw_1x1 / 2) .* exp(+1j*pi*F_tau.^2/Kr);%滤波器
%  匹配滤波
S1_ftau_eta = fftshift(fft(fftshift(S_echo, 2), Nrg, 2), 2);
S1_ftau_eta = S1_ftau_eta .* Hf;
S1_tau_eta = fftshift(ifft(fftshift(S1_ftau_eta, 2), Nrg, 2), 2);
%% 方位向傅里叶变换
S2_tau_feta = fftshift(fft(fftshift(S1_tau_eta, 1), Naz, 1), 1);

%% 距离徙动校正RCMC:采用相位补偿法

%虽然Ka是随着R0变化的，但是在相位补偿时需要假设R0是不变的
delta_R = (((lambda * F_eta).^2) .* R_eta_c) ./ (8 * (Vr^2)); %距离徙动表达式
G_rcmc = exp((+4j * pi .* F_tau .* delta_R)./c); %补偿相位
S3_ftau_feta = fftshift(fft(fftshift(S2_tau_feta, 2), Nrg, 2), 2); %在方位向傅里叶变换的基础上进行距离向傅里叶变换

S3_ftau_feta = S3_ftau_feta .* G_rcmc; %与补偿相位相乘
S3_tau_feta_RCMC = fftshift(ifft(fftshift(S3_ftau_feta, 2), Nrg, 2), 2); %距离向傅里叶逆变换

%距离徙动校正结束

%% 方位压缩
%  随着距离向时间变化的最近斜距；c/2是因为距离向上一个时间包含两次电磁波来回
R0_tau_r = (t_tau * c / 2) * cos(theta_a);
Ext_R0_tau_r = repmat(R0_tau_r, Naz, 1); %扩展R0，用于计算变量Ka

%  根据变化的R0计算出相应的Ka矩阵(距离向变化，方位向不变)
Ka = 2 * Vr^2 * cos(theta_a)^2 ./ (lambda * Ext_R0_tau_r);
%  方位向匹配滤波器
Haz = exp(-1j*pi*F_eta.^2./Ka);
Offset = exp(-1j*2*pi*F_eta.*t_eta_c);%偏移滤波器，将原点搬移到Naz/2的位置，校准坐标
%  匹配滤波
S4_tau_feta = S3_tau_feta_RCMC .* Haz.*Offset;
S4_tau_eta = fftshift(ifft(fftshift(S4_tau_feta, 1), Naz, 1), 1);


%% 回波图可视化

figure('name', "回波可视化")
subplot(2, 2, 1);
imagesc(real(S_echo));
title('原始仿真信号实部');
xlabel('距离向时间 \tau');
ylabel('方位向时间 \eta');

subplot(2, 2, 2);
imagesc(imag(S_echo));
title('原始仿真信号虚部');
xlabel('距离向时间 \tau');
ylabel('方位向时间 \eta');


subplot(2, 2, 3);
imagesc(abs(S_echo));
title('原始仿真信号幅度');
xlabel('距离向时间 \tau');
ylabel('方位向时间 \eta');

subplot(2, 2, 4);
imagesc(angle(S_echo));
title('原始仿真信号相位');
xlabel('距离向时间 \tau');
ylabel('方位向时间 \eta');

%% 距离压缩可视化

%距离频谱可视化
figure('name', "回波频谱可视化")
subplot(2, 2, 1);
imagesc(real(S1_ftau_eta));
title('回波距离向实部');
xlabel('距离向频谱点数f_\tau');
ylabel('方位向时间 \eta');

subplot(2, 2, 2);
imagesc(imag(S1_ftau_eta));
title('回波距离向虚部');
xlabel('距离向频谱点数f_\tau');
ylabel('方位向时间 \eta');

subplot(2, 2, 3);
imagesc(abs(S1_ftau_eta));
title('回波距离向频谱幅度');
xlabel('距离向频谱点数f_\tau');
ylabel('方位向时间 \eta');

subplot(2, 2, 4);
imagesc(angle(S1_ftau_eta));
title('回波距离向相位');
xlabel('距离向频谱点数f_\tau');
ylabel('方位向时间 \eta');

%距离压缩结果
figure('name', "距离压缩时域结果")
subplot(1, 2, 1);
imagesc(real(S1_tau_eta));
title('实部');
xlabel('距离向时间 \tau');
ylabel('方位向时间 \eta');

subplot(1, 2, 2);
imagesc(abs(S1_tau_eta));
title('幅度');
xlabel('距离向时间 \tau');
ylabel('方位向时间 \eta');
%% 方位向傅里叶变换可视化

figure('name', "方位向傅里叶变换结果")
subplot(1, 2, 1);
imagesc(real(S2_tau_feta));
title('实部');
xlabel('距离向时间\tau');
ylabel('方位向频谱点数f_\eta');

subplot(1, 2, 2);
imagesc(abs(S2_tau_feta));
title('幅度');
xlabel('距离向时间\tau');
ylabel('方位向频谱点数f_\eta');

%% 距离徙动校正可视化

figure('name', "距离徙动校正结果")
subplot(1, 2, 1);
imagesc(real(S3_tau_feta_RCMC));
title('实部');
xlabel('距离向时间\tau');
ylabel('方位向频率点数 f_\eta');

subplot(1, 2, 2);
imagesc(abs(S3_tau_feta_RCMC));
title('幅度');
xlabel('距离向时间\tau');
ylabel('方位向频率点数 f_\eta');

%% 回波成像
figure('name', "点目标成像结果")
subplot(1, 2, 1);
imagesc(abs(S4_tau_eta));
title('幅度');
xlabel('距离向时间\tau');
ylabel('方位向时间 \eta');


