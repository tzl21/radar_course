close all;clear all;clc

T = 7.24e-6;                % 信号持续时间
B = 14.2e6;                 % 信号带宽
K = B/T;                    % 调频率
ratio = 1.25;               % 过采样率
Fs = ratio*B;               % 采样频率
dt = 1/Fs;                  % 采样间隔
N = ceil(T/dt);             % 采样点数
t = ((0:N-1)-N/2)/N*T;      % 时间轴
f = ((0:N-1)-N/2)/N*Fs;     % 频率轴

st = exp(1i*pi*K*(t).^2);     % 生成信号
Sf = fft(st);               % FFT
Hf = exp(1i*pi*f.^2/K);     % 频域匹配滤波器

Out = Sf.*Hf;               % 频域匹配滤波


figure,set(gcf,'Color','w');
subplot(1,2,1),plot(f,real(Out));axis tight
title('频谱实部'),xlabel('频率'),ylabel('幅度');
subplot(1,2,2),plot(f,imag(Out));axis tight
title('频谱虚部'),xlabel('频率'),ylabel('幅度');
sgtitle('图3.9 匹配滤波后的信号频谱');