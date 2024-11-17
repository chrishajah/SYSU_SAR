clear;
%%参数设置
bw = 5e6;       %线性调频带宽10Mhz
tw = 20e-6;     %脉宽20us
a = 20;         %过采样率
fs = bw * a;    %采样频率
Ts = 1/fs;      %采样周期
k = bw/tw;      %k

%LFM波形生成
n = -tw/2:Ts:tw/2;
lfm_waveform = exp(1i*pi*k*n.^2);
%定义匹配滤波器
Gf0 = fftshift(fft(lfm_waveform));
Hf = conj(Gf0);
%生成脉冲压缩波形
cs_lfm_waveform = fftshift(ifft(Gf0 .* Hf));
figure(1);
subplot(211);
plot(n,real(lfm_waveform));
ylabel("Amptitude");
xlabel("time/us");
title("原始LFM信号");
xlim([-10e-6,10e-6]);
ylim([-1.2,1.2]);
subplot(212);
plot(n,real(cs_lfm_waveform));
ylabel("Amptitude");
xlabel("time/us");
title("匹配滤波器进行脉冲压缩后的信号");
xlim([-10e-6,10e-6]);
figure(2);
db_cs_lfm_waveform = db(cs_lfm_waveform/max(cs_lfm_waveform));
plot(n,db_cs_lfm_waveform);
xlim([-5e-6,5e-6]);
ylabel("Amptitude/dB");
xlabel("time/us");
title("脉冲压缩信号取dB并放大后的图像");
disp("矩形窗下的脉冲压缩信号：");
disp("IRW = "+num2str(calculate_irw(db_cs_lfm_waveform)));
disp("PSLR = "+num2str(calculate_PSLR(db_cs_lfm_waveform)));
disp("ISLR = "+num2str(calculate_ISLR(cs_lfm_waveform)));
text(-4e-6,-10,"IRW = "+num2str(calculate_irw(db_cs_lfm_waveform))+"us");
text(-4e-6,-15,"PSLR = "+num2str(calculate_PSLR(db_cs_lfm_waveform))+"dB");
text(-4e-6,-20,"ISLR = "+num2str(calculate_ISLR(cs_lfm_waveform))+"dB");

hanning_window = hann(2001).';
Gf1 = fft(hanning_window .* lfm_waveform);
Hf_han = conj(Gf1);
han_cs_lfm_waveform = fftshift(ifft(Hf_han .* Gf1));
db_han_cs_lfm_waveform = db(han_cs_lfm_waveform/max(han_cs_lfm_waveform));
figure(3);
plot(n,db_han_cs_lfm_waveform,'Color','b','DisplayName','Hanning Window');
hold on;
plot(n,db_cs_lfm_waveform,'Color','g','DisplayName','No Window');
legend;


xlim([-5e-6,5e-6]);
ylabel("Amptitude/dB");
xlabel("time/us");
disp("汉宁窗下的脉冲压缩信号：");
disp("IRW = "+num2str(calculate_irw(db_han_cs_lfm_waveform)));
disp("PSLR = "+num2str(calculate_PSLR(db_han_cs_lfm_waveform)));
disp("ISLR = "+num2str(calculate_ISLR(han_cs_lfm_waveform)));
title("加汉宁窗的脉冲压缩信号(in dB)");
text(-4e-6,-10,"Hanning:IRW = "+num2str(calculate_irw(db_han_cs_lfm_waveform))+"us");
text(-4e-6,-15,"PSLR = "+num2str(calculate_PSLR(db_han_cs_lfm_waveform))+"dB");
text(-4e-6,-20,"ISLR = "+num2str(calculate_ISLR(han_cs_lfm_waveform))+"dB");

hanning_window = hamming(2001).';
Gf1 = fft(hanning_window .* lfm_waveform);
Hf_han = conj(Gf1);
han_cs_lfm_waveform = fftshift(ifft(Hf_han .* Gf1));
db_han_cs_lfm_waveform = db(han_cs_lfm_waveform/max(han_cs_lfm_waveform));
plot(n,db_han_cs_lfm_waveform,'DisplayName','Hamming Window');
hold on;
legend;

xlim([-5e-6,5e-6]);
ylabel("Amptitude/dB");
xlabel("time/us");
disp("汉明窗下的脉冲压缩信号：");
disp("IRW = "+num2str(calculate_irw(db_han_cs_lfm_waveform)));
disp("PSLR = "+num2str(calculate_PSLR(db_han_cs_lfm_waveform)));
disp("ISLR = "+num2str(calculate_ISLR(han_cs_lfm_waveform)));
title("加汉明窗的脉冲压缩信号(in dB)");
text(-4e-6,-25,"Hamming :IRW = "+num2str(calculate_irw(db_han_cs_lfm_waveform))+"us");
text(-4e-6,-30,"PSLR = "+num2str(calculate_PSLR(db_han_cs_lfm_waveform))+"dB");
text(-4e-6,-35,"ISLR = "+num2str(calculate_ISLR(han_cs_lfm_waveform))+"dB");


hanning_window = blackman(2001).';
Gf1 = fft(hanning_window .* lfm_waveform);
Hf_han = conj(Gf1);
han_cs_lfm_waveform = fftshift(ifft(Hf_han .* Gf1));
db_han_cs_lfm_waveform = db(han_cs_lfm_waveform/max(han_cs_lfm_waveform));

plot(n,db_han_cs_lfm_waveform,'DisplayName','Blackman Window');
hold on;
legend;

xlim([-5e-6,5e-6]);
ylabel("Amptitude/dB");
xlabel("time/us");
disp("blackman窗下的脉冲压缩信号：");
disp("IRW = "+num2str(calculate_irw(db_han_cs_lfm_waveform)));
disp("PSLR = "+num2str(calculate_PSLR(db_han_cs_lfm_waveform)));
disp("ISLR = "+num2str(calculate_ISLR(han_cs_lfm_waveform)));
title("加blackman窗的脉冲压缩信号(in dB)");
text(-4e-6,-40,"Blackman:IRW = "+num2str(calculate_irw(db_han_cs_lfm_waveform))+"us");
text(-4e-6,-45,"PSLR = "+num2str(calculate_PSLR(db_han_cs_lfm_waveform))+"dB");
text(-4e-6,-50,"ISLR = "+num2str(calculate_ISLR(han_cs_lfm_waveform))+"dB");


hanning_window = kaiser(2001).';
Gf1 = fft(hanning_window .* lfm_waveform);
Hf_han = conj(Gf1);
han_cs_lfm_waveform = fftshift(ifft(Hf_han .* Gf1));
db_han_cs_lfm_waveform = db(han_cs_lfm_waveform/max(han_cs_lfm_waveform));

plot(n,db_han_cs_lfm_waveform,'DisplayName','Kaiser Window');
hold on;
legend;

xlim([-5e-6,5e-6]);
ylabel("Amptitude/dB");
xlabel("time/us");
disp("Kaiser窗下的脉冲压缩信号：");
disp("IRW = "+num2str(calculate_irw(db_han_cs_lfm_waveform)));
disp("PSLR = "+num2str(calculate_PSLR(db_han_cs_lfm_waveform)));
disp("ISLR = "+num2str(calculate_ISLR(han_cs_lfm_waveform)));
title("加窗的脉冲压缩信号(in dB)");
text(-4e-6,-55,"Kaiser:IRW = "+num2str(calculate_irw(db_han_cs_lfm_waveform))+"us");
text(-4e-6,-60,"PSLR = "+num2str(calculate_PSLR(db_han_cs_lfm_waveform))+"dB");
text(-4e-6,-65,"ISLR = "+num2str(calculate_ISLR(han_cs_lfm_waveform))+"dB");



%% 函数实现代码
%% IRW 冲激响应的3dB主瓣宽度
function [irw] = calculate_irw(Af)
    % 找到Af的最大位置
    [~,locmax] = max(Af);
    % 找到locmax左边最接近-3dB的位置
    [~,locleft] = min(abs(Af(1:locmax)+3));
    % 找到locmax右边最接近-3dB的位置
    [~,locright] = min(abs(Af(locmax:end)+3));
    locright = locright + locmax - 1;
    % 得到3dB波束宽度
    irw = locright-locleft;
end
%% PSLR函数 峰值旁瓣比，最大旁瓣与峰值的高度比
function [PSLR] = calculate_PSLR(Af)
    % 找到所有的pesks
    peaks = findpeaks(Af);
    % 对peaks进行降序排列
    peaks = sort(peaks,'descend');
    % 得到第一旁瓣
    PSLR = peaks(2);
end
%% ISLR 实现计算ISLR
function [Ratio] = calculate_ISLR(s)
[M,N] = size(s);
if N==1
    s = s.';
else 
    s = s;
end
L = length(s);
x = abs(s).^2;
ss = x;
[peak,position] = max(x);
 for loop = position:-1:1
     if x(loop)>x(loop-1)
        x(loop) = 0;
     else
        num = loop;
        break; 
     end
end
for loop = position+1:L
    if x(loop)>x(loop+1)
       x(loop) = 0;
    else
        Num = loop;
       break; 
    end
end
y = zeros(1,L);
y(1,num:1:Num) = ss(1,num:1:Num);
Ratio = 10*log10(abs(sum(ss,2)-sum(y,2))/sum(y,2));
end