close all;clear all;clc

T = 42e-6;                % 信号持续时间
B = 17.2e6;                 % 信号带宽
K = B/T;                    % 调频率
ratio = 1.07;               % 过采样率
Fs = ratio*B;               % 采样频率
dt = 1/Fs;                  % 采样间隔
N = ceil(T/dt);             % 采样点数
t = ((0:N-1)-N/2)/N*T;      % 时间轴
f = ((0:N-1)-N/2)/N*Fs;     % 频率轴

st = exp(1i*pi*K*t.^2);     % 生成信号
Sf = fft(st);               % FFT
window = kaiser(N,2.5)';    % 时域窗
Window = fftshift(window);  % 频域窗
Out =exp(1j*pi*f.^2/K);
Hf = Window.*Out;           % 方式3的匹配滤波器：直接在频域生成匹配滤波器


figure,set(gcf,'Color','w');
subplot(3,1,1),plot(abs(Hf));axis tight
title('频域匹配滤波器的幅度'),ylabel('幅度');
subplot(3,1,2),plot(1:1:387,f(1:387)*1e-6+10,'b',388:1:773,f(388:773)*1e-6-10,'b');axis tight
title('频域匹配滤波器的瞬时频率'),ylabel('频率（MHz）');
subplot(3,1,3),plot(fftshift(pi*f.^2/K));axis tight
title('频域匹配滤波器的相位'),ylabel('相位（弧度）'),xlabel('频率（FFT点数）');
sgtitle('图3.12 方式3生成的匹配滤波器');