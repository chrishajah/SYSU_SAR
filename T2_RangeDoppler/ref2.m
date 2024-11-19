%% 合成孔径雷达 SAR回波信号仿真
%% （1）输入参数：
% 场景参数：场景中心距离雷达为10Km，方位×距离——120m×520m，雷达高度为0。
% 目标参数：25个RCS等于1的理想点目标，均匀分布在以场景中心为中心的方位-50m～+50m、距离-250m～+250m的范围内。
% LFM信号参数：中心频率10.0GHz，脉冲宽度30us。
% SAR的方位斜视角等于0，即正侧视。
% 天线加权形式：矩形加权或sinc加权。
% 分辨率有两种情况：（1）10m×10m；（2）1m×1m。
%% （2）针对两种分辨率，完成回波仿真：
% 设计发射LFM信号带宽，保证图像的分辨率要求，同时又不过分增加仿真的数据量。
% 设计天线尺寸，保证距离向覆盖场景信号，图像的分辨率要求，同时又不过分增加仿真的数据量。
% 设计雷达接收机波门起始时刻和中止时刻，保证场景中的每一个目标的回波都会被完全接收，同时又不过分增加仿真的数据量。
% 设计两个方向的过采样因子（最好一致），使得距离和方位向的采样模糊可以忽略，同时又不过分增加仿真的数据量。
% 按正交解调接收方式仿真SAR回波信号，观察信号的时、频域波形。
%% （3）完成对回波信号的分析，理解方位向信号的特性：
% 分析回波信号在方位向的时频分布，从信号角度理解方位向高分辨率原理。
% 分析回波信号在二维时域、二维频谱、距离多普勒域的特性。
% 分析距离压缩后的回波数据，理解距离迁移（RCMC）的特性。
%% （4）对每个距离向脉冲完成脉冲压缩
clear;clc;close all;
%% 定义参数
% --------------------------------------------------------------------
R_eta_c = 10e3;                % 景中心距10km
H = 0;                      % 雷达高度为0
V_r = 340;                  % 雷达有效速度，因为雷达高度为0，所以 波束速度Vg=平台速度Vs
T_r = 30e-6;                % 发射脉冲时宽
f_0 = 10e9;                 % 雷达工作频率
rho = 10;                   % 分辨率有两种情况：(1)10m×10m  (2)1m×1m
rho = 1;
C = 3e8;                    % 光速
%% 设计发射LFM信号带宽，保证图像的分辨率要求，同时又不过分增加仿真的数据量。
BW_range = C/2/rho;         % 距离向带宽 BW_range = Kr*Tr;       
K_r = BW_range/T_r;         % 距离调频率
%% 设计两个方向的过采样因子（最好一致），使得距离和方位向的采样模糊可以忽略，同时又不过分增加仿真的数据量。
alpha_s_r = 1.2;             % 距离向过采样因子
alpha_s_a = 1.2;             % 方位向过采样因子
% 方位向过采样率通常为1.1~1.4.如果PRF太低，由混叠引起的方位模糊会很严重
% 方位向过采样率通常大于距离向过采样率，因为方位频谱比距离频谱衰落得慢
if rho==1
    N_azimuth = 256;          	% 距离线数（即数据矩阵，行数）   
    N_range = 6500;             % 距离线采样点数（即数据矩阵，列数）
    NFFT_azimuth = 256;         % 方位向FFT长度    
    NFFT_range = 8192;          % 距离向FFT长度
else
    N_azimuth = 32;          	% 距离线数（即数据矩阵，行数）   
    N_range = 650;              % 距离线采样点数（即数据矩阵，列数）
    NFFT_azimuth = 32;          % 方位向FFT长度    
    NFFT_range = 1024;          % 距离向FFT长度
end
% NFFT_r = Nrg;             % 距离向FFT长度
% NFFT_a = Naz;             % 方位向FFT长度
theta_r_c = (0*pi)/180;	    % 波束斜视角，0 度，转换为弧度
R_0 = sqrt((R_eta_c*cos(theta_r_c))^2+H^2);  % 假设的雷达高度为0，高度又为H=0，雷达到目标的最短距离记为R0
lambda = C/f_0;             % 波长
f_eta_c = 2*V_r*sin(theta_r_c)/lambda;% 多普勒中心频率，根据公式（4.33）计算
F_r = alpha_s_r*BW_range;     % 距离采样率
% alpha_sr = Fr / BW_range;   % 距离向过采样因子  Fr/Kr/Tr
% alpha_sa = Fa / BW_dop;     % 方位向过采样因子  Fa/BW_dop
Nr = T_r*F_r;                % 线性调频信号采样点数
K_a = 2*V_r^2*(cos(theta_r_c))^3./(lambda.* R_0);  	
% 方位向调频率---是方位向频率或多普勒频率的变化率，根据公式（4.38）计算的
%% 设计天线尺寸，保证距离向覆盖场景信号，图像的分辨率要求，同时又不过分增加仿真的数据量。
La_real = 2*rho;   
% if rho==1
%     La_real = 2;                                    % 方位向天线长度
%     % 能达到的方位向分辨率约为天线实际长度的一半
%     % La_real = 0.886*2*V_r*cos(theta_r_c)/BW_dop;   % 根据公式（4.36）
% else
%     La_real = 20;                                   % 方位向天线长度
%     % La_real = 0.886*2*V_r*cos(theta_r_c)/BW_dop;   % 根据公式（4.36）
% end
beta_bw = 0.886*lambda/La_real; % 雷达3dB波束
L_a = beta_bw*R_0;              % 合成孔径长度
%% 多普勒带宽
BW_dop = 0.886*2*V_r*cos(theta_r_c)/La_real;     % 根据公式（4.36）
F_a = alpha_s_a*BW_dop;    % 方位采样率;方位谱衰减得比较慢，过采样率一般取1.1~1.4，以减小方位模糊功率
Mamb = round(f_eta_c/F_a);    % 多普勒模糊
% --------------------------------------------------------------------
%% 设定仿真点目标的位置,以距离向作为x轴正方向,以方位向作为y轴正方向
% 25个RCS等于1的理想点目标，均匀分布在以场景中心为中心的方位-50m～+50m、距离-250m～+250m的范围内
Ptar = [R_0-250,-50;R_0-250,-25;R_0-250,0;R_0-250,25;R_0-250,50;
      R_0-125,-50;R_0-125,-25;R_0-125,0;R_0-125,25;R_0-125,50;
      R_0,-50;R_0,-25;R_0,0;R_0,25;R_0,50;
      R_0+125,-50;R_0+125,-25;R_0+125,0;R_0+125,25;R_0+125,50;
      R_0+250,-50;R_0+250,-25;R_0+250,0;R_0+250,25;R_0+250,50;];%距离向坐标,方位向坐标
[Ntar,n] = size(Ptar);

eta_c_target = zeros(Ntar,1);% 计算目标各自的波束中心穿越时刻
Target_x_range = zeros(Ntar,1);
Target_y_azimuth = zeros(Ntar,1);
for i=1:Ntar
    Target_x_range(i) = Ptar(i,1);
    Target_y_azimuth(i) = Ptar(i,2);
    eta_c_target(i) = (Target_y_azimuth(i)-Target_x_range(i)*tan(theta_r_c))/V_r;% 目标i的波束中心穿越时刻
end
%% 
% 距离（方位）向时间，频率相关定义
% --------------------------------------------------------------------
% 距离时间也称为"快时间"，距离间隔由距离时间和光速得到
% t_tau = linspace((2*(R_0-250)/C),(2*(R_0+250)/C),6500);
t_tau = 2*R_eta_c/C + ( -N_range/2 : (N_range/2-1) )/F_r;                      % 距离时间轴
fr = ( -NFFT_range/2 : NFFT_range/2-1 )*( F_r/NFFT_range );              % 距离频率轴
% 方位时间也称为"慢时间"，方位间隔由方位时间和比光速慢得多的波束照射区推移速度得到
% 目标方位时间以零多普勒为参考原点。对于多目标情况，应给定一共同的绝对时刻\eta（如数据截获起始时刻\eta=0）
% t_eta = linspace(-(L_a+100)/V_r/2,(L_a+100)/V_r/2-1/F_a,256);
t_eta = ( -N_azimuth/2: N_azimuth/2-1 )/F_a;                               % 方位时间轴
fa = f_eta_c + ( -NFFT_azimuth/2 : NFFT_azimuth/2-1 )*( F_a/NFFT_azimuth );	% 方位频率轴
 
% 生成距离（方位）时间（频率）矩阵
tau_mtx_range = ones(N_azimuth,1)*t_tau;        % 距离时间轴矩阵，大小：Naz*Nrg
tau_mtx_azimuth = t_eta.'*ones(1,N_range);      % 方位时间轴矩阵，大小：Naz*Nrg
%% 
% --------------------------------------------------------------------
% 生成点目标原始数据
% --------------------------------------------------------------------
s_0_tau_eta = zeros(N_azimuth,N_range);      % 用来存放生成的回波数据，所有点目标回波信号之和  
RCS = 1;                                 % 理想点目标，RCS均为1，目标回波幅度都设置为1.
R_eta = zeros(N_azimuth,N_range,Ntar);
w_range = zeros(N_azimuth,N_range);
si_tau_eta = zeros(N_azimuth,N_range,Ntar);
for k = 1:Ntar                          % 生成k个目标的原始回波数据
    R_eta(:,:,k) = sqrt( (Target_x_range(k).*ones(N_azimuth,N_range)).^2 + (V_r.*tau_mtx_azimuth-Target_y_azimuth(k).*ones(N_azimuth,N_range)).^2 );% 目标k的瞬时斜距
    w_range = ((abs(tau_mtx_range-2.*R_eta(:,:,k)./C)) <= ((T_r/2).*ones(N_azimuth,N_range)));     % 距离向包络，即距离窗
    % =====================================================================    
    % 方位向包络，也就是 天线的双程方向图作用因子。
    % 方式1
    % sinc平方型函数，根据公式（4.31）计算    
    theta = atan( V_r.*(tau_mtx_azimuth-eta_c_target(k).*ones(N_azimuth,N_range))/Target_x_range(k) );
    w_azimuth = (sinc(0.886.*theta./beta_bw)).^2;    
    % 用每个目标对应的 波束中心穿越时刻，而不是之前参数中的nc。
    
    % 方式2
    % 利用合成孔径长度，直接构造矩形窗（其实这里只是限制数据范围，没有真正加窗）
    % 天线加权形式：矩形加权
%     w_azimuth = (abs(t_eta - eta_c_target(k)) <= (L_a/2)/V_r);    % 行向量
%     w_azimuth = w_azimuth.'*ones(1,N_range);    % 生成Naz*Nrg的矩阵
    % =====================================================================     
    % s_r = RCS.*w_range.*w_azimuth.*cos(2*pi*f_0*(tau_mtx_range-2*R_eta(:,:,k)./C)+pi*K_r*(tau_mtx_range-2*R_eta(:,:,k)./C).^2)+Psi;
    % 上式就是生成的某一个点目标（目标k）的回波信号,Psi表示地表散射过程可能引起的雷达信号相位改变。
    s_0_tau_eta_k = RCS.*w_range.*w_azimuth.*exp(-(1j*4*pi*f_0).*R_eta(:,:,k)./C).*exp((1j*pi*K_r).*(tau_mtx_range-2.*R_eta(:,:,k)./C).^2);
    % 上式是由式（4.39）给出的解调后单个点目标的基带信号复数表达式
    % 经过25次循环，生成25个点目标的回波信号，相加即可
    si_tau_eta(:,:,k) = s_0_tau_eta_k;
    s_0_tau_eta = s_0_tau_eta + s_0_tau_eta_k;  % 所有点目标回波信号之和（已正交解调）   
end
% s_0_tau_eta 就是点目标回波信号原始数据
%%  作图
% 按正交解调接收方式仿真SAR回波信号，观察信号的时、频域波形。
figure(1); % 时域
subplot(2,2,1);
imagesc(real(s_0_tau_eta));
title('实部');xlabel('距离时域（采样点）');ylabel('方位时域（采样点）');
subplot(2,2,2);
imagesc(imag(s_0_tau_eta));
title('虚部');xlabel('距离时域（采样点）');ylabel('方位时域（采样点）');
subplot(2,2,3);
imagesc(abs(s_0_tau_eta));
title('幅度');xlabel('距离时域（采样点）');ylabel('方位时域（采样点）');
subplot(2,2,4);
imagesc(angle(s_0_tau_eta));
title('相位');xlabel('距离时域（采样点）');ylabel('方位时域（采样点）');
%% 作图---观察距离徙动现象
figure(2);
subplot(2,1,1);
imagesc(real(si_tau_eta(:,:,13))); xlim([2500 4000])
title('实部');xlabel('距离时域（采样点）');ylabel('方位时域（采样点）');
subplot(2,1,2);
imagesc(imag(si_tau_eta(:,:,13))); xlim([2500 4000])
title('虚部');xlabel('距离时域（采样点）');ylabel('方位时域（采样点）');
figure(3);
subplot(2,1,1);
imagesc(abs(si_tau_eta(:,:,13)));xlim([2500 4000])
title('幅度');xlabel('距离时域（采样点）');ylabel('方位时域（采样点）');
subplot(2,1,2);
imagesc(angle(si_tau_eta(:,:,13)));xlim([2500 4000])
title('相位');xlabel('距离时域（采样点）');ylabel('方位时域（采样点）');
%% 按正交解调接收方式仿真SAR回波信号，观察信号的时、频域波形。
figure(4);% 频域
subplot(2,2,1);
imagesc(abs(fft(s_0_tau_eta))); % 距离时域-方位频域
xlabel('距离向(采样点)'),ylabel('方位频率(采样点)'),title('RD频谱幅度');
subplot(2,2,2);
imagesc(angle(fft(s_0_tau_eta)));
xlabel('距离向(采样点)'),title('RD频谱相位');
subplot(2,2,3);
imagesc(abs(fft2(s_0_tau_eta))); % 距离频域-方位频域
xlabel('距离频率(采样点)'),ylabel('方位频率(采样点)'),title('二维频谱幅度');
subplot(2,2,4);
imagesc(angle(fft2(s_0_tau_eta)));% fft2 二维快速傅里叶变换
xlabel('距离频率(采样点)'),title('二维频谱相位');
%% 距离压缩
% --------------------------------------------------------------------
S_range = fft(s_0_tau_eta,NFFT_range,2);     % 进行距离向傅里叶变换，零频在两端。
% 作图
% 图——距离频域，方位时域，频谱（未距离压缩）
figure(5);
subplot(1,2,1);
imagesc(real(S_range)); % xlim([6000 6200]);
title('实部');xlabel('距离频域（采样点）');ylabel('方位时域（采样点）'); 
subplot(1,2,2);
imagesc(abs(S_range)); % xlim([6000 6200]);
title('幅度');xlabel('距离频域（采样点）');ylabel('方位时域（采样点）');
% 生成距离向匹配滤波器
% 采用方式2---通过距离多普勒域插值实现随距离的变化的距离徙动校正
% 时域复制脉冲，末端补零，fft，再取复共轭。
t_ref = ( -Nr/2 : (Nr/2-1) )/F_r;               % 用来生成距离MF的距离时间轴
t_ref_mtx = ones(N_azimuth,1)*t_ref;            % 矩阵形式
w_ref = kaiser(Nr,2.5);                         % 距离向，构建Kaiser窗，此为列向量。
w_ref = ones(N_azimuth,1)*(w_ref.');            % 构成矩阵形式，每一行都相同的加窗。
s_ref = exp((1j*pi*K_r).*((t_ref_mtx).^2));     % 复制（发射）脉冲，未加窗。
% s_ref = w_ref.*exp((1j*pi*Kr).*((t_ref_mtx).^2)); % 复制（发射）脉冲，加了窗。
s_ref = [s_ref,zeros(N_azimuth,N_range-Nr)];    % 对复制脉冲，后端补零。
S_ref = fft(s_ref,NFFT_range,2);                % 复制脉冲的距离傅里叶变换，零频在两端。
H_f_tau = conj(S_ref);                          % 距离向匹配滤波器，零频在两端。
S_range_c = S_range.*H_f_tau;                   % 乘以匹配滤波器，零频在两端。    
s_rc_tau_eta = ifft(S_range_c,[],2);            % 完成距离压缩，回到二维时域。s_rc_tau_eta即为距离压缩输出
% s_rc的长度为：Naz*Nrg。未去除弃置区。对s_rc进行去除弃置区的操作
% 弃置区长度为：2*（Nr-1）
% 我们截取的长度：（Nrg-Nr+1），记为 N_rg。
N_rg = N_range-Nr+1;                        % 完全卷积的长度
% s_rc_c = zeros(Naz,N_rg);               % 用来存放去除弃置区后的数据
s_rc_c = s_rc_tau_eta(:,1:N_rg);                % 取前 N_rg列。
% ====================================================
% 作图——距离频域，方位时域，频谱（已距离压缩）
figure(6);
subplot(1,2,1);
imagesc(real(S_range_c));% xlim([6000 6200]);
title('实部');xlabel('距离向（采样点）');ylabel('方位向（采样点）');     
subplot(1,2,2);
imagesc(abs(S_range_c));% xlim([6000 6200]);
title('幅度');xlabel('距离向（采样点）');ylabel('方位向（采样点）');
% 作图——二维时域（完成距离压缩）
figure(7);
subplot(1,2,1);
imagesc(real(s_rc_c));  %　这及其以下，都直接使用去除弃置区后的结果
title('实部');xlabel('距离时域（采样点）');ylabel('方位时域（采样点）');
subplot(1,2,2);
imagesc(abs(s_rc_c));% xlim([395 410]);
title('幅度');xlabel('距离时域（采样点）');ylabel('方位时域（采样点）');
%% 变换到距离多普勒域，进行距离徙动校正
s_rc_c = s_rc_c.*exp(-1j*2*pi*f_eta_c.*(t_eta.'*ones(1,N_rg)));    % 数据搬移
S_rd = fft(s_rc_c,NFFT_azimuth,1);            % 方位向傅里叶变换，到距离多普勒域
% ====================================================================
% 设置方位频率轴
fa = f_eta_c + fftshift(-NFFT_azimuth/2:NFFT_azimuth/2-1)/NFFT_azimuth*F_a;    % 方位频率轴如此设置。
% =====================================================================
% 下面这个是改进的:每一个最近斜距（R0）都随着距离门的不同而改变。
tr_RCMC = 2*Target_x_range(13)/C + ( -N_rg/2 : (N_rg/2-1) )/F_r;   % 在新的距离线长度下的时间轴。
R0_RCMC = (C/2).*tr_RCMC;       % 随距离线变化的R0，记为R0_RCMC，用来计算RCM和Ka。
delta_Rrd_fn = lambda^2.*((fa.').^2)*(R0_RCMC)/(8*V_r^2);
num_range = C/(2*F_r);   % 一个距离采样单元，对应的长度
delta_Rrd_fn_num = delta_Rrd_fn./num_range; % 每一个方位向频率，其RCM对应的距离采样单元数
R = 8;  % sinc插值核长度
S_rd_rcmc = zeros(NFFT_azimuth,N_rg); % 用来存放RCMC后的值
for p = 1 : NFFT_azimuth
    for q = 1 : N_rg   % 此时距离向的长度是 (Nrg-Nr+1)=N_rg        
        delta_Rrd_fn_p = delta_Rrd_fn_num(p,q);
        Rrd_fn_p = q + delta_Rrd_fn_p;
        Rrd_fn_p_zheng = ceil(Rrd_fn_p);        % ceil，向上取整。
        ii = ( Rrd_fn_p-(Rrd_fn_p_zheng-R/2):-1:Rrd_fn_p-(Rrd_fn_p_zheng+R/2-1)  );        
        rcmc_sinc = sinc(ii);
        rcmc_sinc = rcmc_sinc/sum(rcmc_sinc);   % 插值核的归一化
        % ii 是sinc插值过程的变量;
        % g(x)=sum(h(ii)*g_d(x-ii)) = sum(h(ii)*g_d(ll));     
        % 由于S_rd只有整数点取值，且范围有限。因此插值中要考虑它的取值溢出边界问题。
        % 这里我采取循环移位的思想，用来解决取值溢出问题。
        if (Rrd_fn_p_zheng-R/2) > N_rg    % 全右溢
            ll = (Rrd_fn_p_zheng-R/2-N_rg:1:Rrd_fn_p_zheng+R/2-1-N_rg);
        else
            if (Rrd_fn_p_zheng+R/2-1) > N_rg    % 部分右溢
                ll_1 = (Rrd_fn_p_zheng-R/2:1:N_rg);
                ll_2 = (1:1:Rrd_fn_p_zheng+R/2-1-N_rg);
                ll = [ll_1,ll_2];
            else
                if (Rrd_fn_p_zheng+R/2-1) < 1    % 全左溢（不可能发生，但还是要考虑）
                    ll = (Rrd_fn_p_zheng-R/2+N_rg:1:Rrd_fn_p_zheng+R/2-1+N_rg);
                else
                    if (Rrd_fn_p_zheng-R/2) < 1       % 部分左溢
                        ll_1 = (Rrd_fn_p_zheng-R/2+N_rg:1:N_rg);
                        ll_2 = (1:1:Rrd_fn_p_zheng+R/2-1);
                        ll = [ll_1,ll_2];
                    else
                        ll = (Rrd_fn_p_zheng-R/2:1:Rrd_fn_p_zheng+R/2-1);
                    end                    
                end
            end
        end   
        rcmc_S_rd = S_rd(p,ll);
        S_rd_rcmc(p,q) = sum( rcmc_sinc.*rcmc_S_rd );
    end
end
% S_rd_rcmc 就是RCMC后的距离多普勒域频谱。
% 绘图——距离多普勒域（未RCMC）
figure(8);
subplot(1,2,1);
imagesc(real(S_rd));
title('实部');xlabel('距离时域（采样点）');ylabel('方位频域（采样点）');      
subplot(1,2,2);
imagesc(abs(S_rd)); %xlim([395 410]);
title('幅度');xlabel('距离时域（采样点）');ylabel('方位频域（采样点）');
%% 绘图——距离多普勒域，RCMC后的结果.RCMC将距离徙动曲线拉直到与方位频率轴平行的方向。
figure(9);
subplot(1,2,1);
imagesc(real(S_rd_rcmc));
title('实部');xlabel('距离时域（采样点）');ylabel('方位频域（采样点）');    
subplot(1,2,2);
imagesc(abs(S_rd_rcmc));
title('幅度');xlabel('距离时域（采样点）');ylabel('方位频域（采样点）');
%% --------------------------------------------------------------------
% 方位压缩---通过每一距离门上的频域匹配滤波实现方位压缩
fa_azimuth_MF = fa;         % 方位频率轴，采用和RCMC中所用的频率轴相同。
K_a = 2*V_r^2*(cos(theta_r_c))^3./(lambda.* R0_RCMC);  	% 方位向调频率，是随最近斜距R0变化的。
Ka_1 = 1./K_a;                                       % 为了方便计算，取方位向调频率的倒数。
Haz = exp( -1j*pi.*(((fa_azimuth_MF).').^2*Ka_1) );	% 方位向匹配滤波器
S_rd_c = S_rd_rcmc.*Haz;            % 乘以匹配滤波器
% 最后通过方位IFFT将数据变换回时域，得到压缩后的复图像。（如果可以还能增加幅度检测及多视叠加）
s_ac = ifft(S_rd_c,[],1);       	% 完成方位压缩，变到图像域。结束。
% 绘图 成像结果
figure(10);
imagesc(abs(s_ac));
title('点目标成像');xlabel('距离时域（采样点）');ylabel('方位时域（采样点）');     
%==========================================================================================================================