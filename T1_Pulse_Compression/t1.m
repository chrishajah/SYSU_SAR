%% 
clear;

%基本参数

bw = 30e6;      %线性调频带宽30Mhz
tw = 30e-6;     %脉宽30us
fs = bw;        %采样频率
Ts = 1/fs;      %采样周期
k = bw/tw;      %k

n = -tw/2:Ts:tw/2;

t = linspace(-30, 30, 901);
lfm_waveform = exp(1i*pi*k*n.^2);
phi = pi * k * n.^2;
ft = k*n;

figure(1);
subplot(221);
plot(t,real(lfm_waveform));
ylabel("Real part of LFM signal");
xlabel("time/us");
xlim([-30,30]);
ylim([-1.2,1.2]);
subplot(222);
plot(t,imag(lfm_waveform));
ylabel("Imaginary part of LFM signal");
xlabel("time/us");
xlim([-30,30]);
ylim([-1.2,1.2]);
subplot(223);
plot(t,phi);
ylabel("信号相位 phi");
xlabel("time/us");
xlim([-30,30]);

subplot(224);
plot(t,ft);
ylabel("瞬时频率 f");
xlabel("time/us");
xlim([-30,30]);

%% 


freq = linspace(-30, 30, 2048);

sample1_fft = fft(lfm_waveform,2048);
sample1_fft = fftshift(sample1_fft);
figure(2);
subplot(222);
plot(freq,abs(sample1_fft));
ylabel("过采样率=1下的信号频谱");
xlabel("Freq/MHz");


fs = bw*1.2;            %采样频率，定义为1.2倍过采样
Ts = 1/fs;              %采样周期
n = -tw/2:Ts:tw/2;
lfm_waveform = exp(1i*pi*k*n.^2);
sample2_fft = fft(lfm_waveform,2048);
sample2_fft = fftshift(sample2_fft);
subplot(223);
plot(freq,abs(sample2_fft));
ylabel("过采样率=1.2下的信号频谱");
xlabel("Freq/MHz");



fs = bw*1.5;            %采样频率，定义为1.5倍过采样
Ts = 1/fs;              %采样周期
n = -tw/2:Ts:tw/2;
lfm_waveform = exp(1i*pi*k*n.^2);
sample3_fft = fft(lfm_waveform,2048);
sample3_fft = fftshift(sample3_fft);
subplot(224);
plot(freq,abs(sample3_fft));
ylabel("过采样率=1.5下的信号频谱");
xlabel("Freq/MHz");


fs = bw*0.8;            %采样频率，定义为0.8倍过采样
Ts = 1/fs;              %采样周期
n = -tw/2:Ts:tw/2;
lfm_waveform = exp(1i*pi*k*n.^2);
sample0_fft = fft(lfm_waveform,2048);
sample0_fft = fftshift(sample0_fft);
subplot(221);
plot(freq,abs(sample0_fft));
ylabel("过采样率=0.8下的信号频谱");
xlabel("Freq/MHz");

%% 


bw = 5e6;           %线性调频带宽5Mhz
tw = 5e-6;          %脉宽5us
a = 2;              %过采样率为2
fs = bw*a;          %采样频率
Ts = 1/fs;          %采样周期
k = bw/tw;          %k
n = -tw/2:Ts:tw/2;
freq = linspace(-fs/2,fs/2,2048);


lfm_waveform = exp(-1i*pi*k*n.^2);
Gf0 = sqrt(tw*fs*a)*exp(-1i*pi*freq.^2/k+1i*pi/4).*(abs(freq)<=(k*tw/2));
Gf1 = fftshift(fft(lfm_waveform,2048));
figure(3);
subplot(521);
plot(freq,abs(Gf1),'DisplayName','FFT');
hold on;
plot(freq,abs(Gf0),'color','r','LineStyle','--','DisplayName','POSP');
grid on;
legend;
xlabel("Freq/Hz");
ylabel("Amptitude");
text(-4.5e6,13,"TBP=25");
subplot(522);
plot(freq,-unwrap(angle(Gf0)) ,'DisplayName','FFT');
hold on;
plot(freq,-((-pi*freq.^2/k)-max(-pi*freq.^2/k))-20,'color','r','LineStyle','--','DisplayName','POSP');
grid on;
xlim([-2e6,2e6]);
legend;
xlabel("Freq/Hz");
ylabel("Phase");
text(-4.5e6,13,"TBP=25");

bw = 5e6;           %线性调频带宽5Mhz
tw = 10e-6;          %脉宽10us
a = 2;              %过采样率为2
fs = bw*a;          %采样频率
Ts = 1/fs;          %采样周期
k = bw/tw;          %k
n = -tw/2:Ts:tw/2;
freq = linspace(-fs/2,fs/2,2048);


lfm_waveform = exp(-1i*pi*k*n.^2);
Gf0 = sqrt(tw*fs*a)*exp(-1i*pi*freq.^2/k+1i*pi/4).*(abs(freq)<=(k*tw/2));
Gf1 = fftshift(fft(lfm_waveform,2048));
figure(3);
subplot(523);
plot(freq,abs(Gf1),'DisplayName','FFT');
hold on;
plot(freq,abs(Gf0),'color','r','LineStyle','--','DisplayName','POSP');
grid on;
legend;
xlabel("Freq/Hz");
ylabel("Amptitude");
text(-4.5e6,13,"TBP=50");
subplot(524);
plot(freq,-unwrap(angle(Gf0)) ,'DisplayName','FFT');
hold on;
plot(freq,-((-pi*freq.^2/k)-max(-pi*freq.^2/k))-38,'color','r','LineStyle','--','DisplayName','POSP');
grid on;
xlim([-2e6,2e6]);
legend;
xlabel("Freq/Hz");
ylabel("Phase");



bw = 10e6;           %线性调频带宽5Mhz
tw = 10e-6;         %脉宽10us
a = 2;              %过采样率为2
fs = bw*a;          %采样频率
Ts = 1/fs;          %采样周期
k = bw/tw;          %k
n = -tw/2:Ts:tw/2;
freq = linspace(-fs/2,fs/2,2048);


lfm_waveform = exp(-1i*pi*k*n.^2);
Gf0 = sqrt(tw*fs*a)*exp(-1i*pi*freq.^2/k+1i*pi/4).*(abs(freq)<=(k*tw/2));
Gf1 = fftshift(fft(lfm_waveform,2048));
figure(3);
subplot(525);
plot(freq,abs(Gf1),'DisplayName','FFT');
hold on;
plot(freq,abs(Gf0),'color','r','LineStyle','--','DisplayName','POSP');
grid on;
legend;
xlabel("Freq/Hz");
ylabel("Amptitude");
text(-4.5e6,13,"TBP=100");
subplot(526);
plot(freq,-unwrap(angle(Gf0)) ,'DisplayName','FFT');
hold on;
plot(freq,-((-pi*freq.^2/k)-max(-pi*freq.^2/k))-76,'color','r','LineStyle','--','DisplayName','POSP');
grid on;
xlim([-5e6,5e6]);
legend;
xlabel("Freq/Hz");
ylabel("Phase");

bw = 20e6;           %线性调频带宽20Mhz
tw = 10e-6;         %脉宽10us
a = 2;              %过采样率为2
fs = bw*a;          %采样频率
Ts = 1/fs;          %采样周期
k = bw/tw;          %k
TBP = k * tw^2;     %TBP=200
n = -tw/2:Ts:tw/2;
freq = linspace(-fs/2,fs/2,2048);


lfm_waveform = exp(-1i*pi*k*n.^2);
Gf0 = sqrt(tw*fs*a)*exp(-1i*pi*freq.^2/k+1i*pi/4).*(abs(freq)<=(k*tw/2));
%驻定相位法
Gf1 = fftshift(fft(lfm_waveform,2048));
%直接FFT
figure(3);
subplot(527);
plot(freq,abs(Gf1),'DisplayName','FFT');
hold on;
plot(freq,abs(Gf0),'color','r','LineStyle','--','DisplayName','POSP');
grid on;
legend;
xlabel("Freq/Hz");
ylabel("Amptitude");
text(-9e6,13,"TBP=200");
subplot(528);
plot(freq,-unwrap(angle(Gf0)) ,'DisplayName','FFT');
hold on;
plot(freq,-((-pi*freq.^2/k)-max(-pi*freq.^2/k))-157,'color','r','LineStyle','--','DisplayName','POSP');
grid on;
xlim([-10e6,10e6]);
legend;
xlabel("Freq/Hz");
ylabel("Phase");

bw = 20e6;           %线性调频带宽20Mhz
tw = 20e-6;         %脉宽10us
a = 2;              %过采样率为2
fs = bw*a;          %采样频率
Ts = 1/fs;          %采样周期
k = bw/tw;          %k
%TBP=400
n = -tw/2:Ts:tw/2;
freq = linspace(-fs/2,fs/2,2048);


lfm_waveform = exp(-1i*pi*k*n.^2);
Gf0 = sqrt(tw*fs*a)*exp(-1i*pi*freq.^2/k+1i*pi/4).*(abs(freq)<=(k*tw/2));
Gf1 = fftshift(fft(lfm_waveform,2048));
figure(3);
subplot(529);
plot(freq,abs(Gf1),'DisplayName','FFT');
hold on;
plot(freq,abs(Gf0),'color','r','LineStyle','--','DisplayName','POSP');
grid on;
legend;
xlabel("Freq/Hz");
ylabel("Amptitude");
text(-9e6,13,"TBP=400");
subplot(5,2,10);
plot(freq,-unwrap(angle(Gf0)) ,'DisplayName','FFT');
hold on;
plot(freq,-((-pi*freq.^2/k)-max(-pi*freq.^2/k))-314,'color','r','LineStyle','--','DisplayName','POSP');
grid on;
xlim([-10e6,10e6]);
legend;
xlabel("Freq/Hz");
ylabel("Phase");