clear;
%%信号波形参数设置
fc = 5e9;       %中心频率5Ghz
bw = 30e6;      %线性调频带宽30Mhz
tw = 30e-6;     %脉宽30us
a = 1.2;        %过采样率
fs = bw * a;    %采样频率
Ts = 1/fs;      %采样周期
k = bw/tw;      %k
c = 3e8;        %光速
%LFM波形生成
n = -tw/2:Ts:tw/2;
lfm_waveform = exp(1i*pi*k*n.^2+1i*pi*fc);

%仿真15、16km回波信号
%利用雷达方程计算目标回波信噪比
t2_snr = radar_eq(5,fc,45,1,290,bw,3,6,16e3);
t1_snr = radar_eq(5,fc,45,1,290,bw,3,6,15e3);
%计算回波延迟
t1_delay = range2time(15e3)-1e-4;
t2_delay = range2time(16e3)-1e-4;
t1_delay_n = round(t1_delay/Ts);
t2_delay_n = round(t2_delay/Ts);


t_rec = linspace(0,0.04,round(0.04e-3/Ts));
%回波接收信号定义
t1_sig = zeros(1,round(0.04e-3/Ts));
for i=t1_delay_n:1:t1_delay_n+length(lfm_waveform)-1
    t1_sig(i+1) = lfm_waveform(i+1 - t1_delay_n);
end

t2_sig = zeros(1,round(0.04e-3/Ts));
for i=t2_delay_n:1:t2_delay_n+length(lfm_waveform)-1
    t2_sig(i+1) = lfm_waveform(i+1 - t2_delay_n);
end

recieve_sig = awgn(t1_sig,t1_snr)+awgn(t2_sig,t2_snr);


figure(1);
plot(t_rec,real(recieve_sig));
ylabel("Amptitude");
xlabel("time/ms");
title("接收原始回波信号");

Gf1 = fft(recieve_sig);
%定义匹配滤波器
Gf0 = fft(real(lfm_waveform),length(Gf1));
Hf = conj(Gf0);
cs_lfm_waveform = ifft(Gf1 .* Hf);
figure(2);
plot(t_rec,abs(cs_lfm_waveform));
ylabel("Amptitude");
xlabel("time/ms");
title("脉冲压缩后的回波信号");

indices = find(abs(cs_lfm_waveform) > 300).*0.04e-3/round(0.04e-3/Ts);
indices = time2range(indices);

%radar  function
function [snr] = radar_eq(pt,freq,g,sigma,te,b,nf,loss,range)
%This is a program of radar eq
c=3.0e+8; %speed of light
lambda =c/freq; %wavelength
p_peak=10*log10(pt);  %convert peak power to dB
lambda_sqdb=10*log10(lambda^2);  %computr wavelength square in dB
sigmadb=10*log10(sigma);%convert sigma to dB
four_pi_cub=10*log10((4*pi)^3);  %(4pi)^3 in dB
k_db=10*log10(1.3e-23);%boltzman's constant in dB
te_db=10*log10(te);  %noisetemp. in dB
b_db=10*log10(b);  %bandwidth in dB
range_pwr4_db=10*log10(range.^4);%vector of target range^4 in dB
%implement Equation(1.56)
num=p_peak+g+g+lambda_sqdb+sigmadb;%分子
den=four_pi_cub+k_db+te_db+b_db+nf+loss+range_pwr4_db;%分母
snr=num-den;
end
