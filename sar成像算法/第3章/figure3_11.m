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
window = kaiser(N,2.5)';    % 时域窗
st_window = st.*window;     % 时域加窗

Hf= conj(fft(st_window,2048));  % FFT频域滤波


figure,set(gcf,'Color','w');
subplot(2,1,1),plot(abs(Hf));axis tight
title('加权后匹配滤波器的幅度谱'),ylabel('幅度');
subplot(2,1,2),plot(1:1:990,-unwrap(angle(Hf(1:990))),'b',1060:1:2048,-unwrap(angle(Hf(1060:2048))),'b');axis tight
title('匹配滤波器的相位谱'),ylabel('相位（弧度）'),xlabel('频率（FFT点数）');
sgtitle('图3.11 方式2生成的匹配滤波器频率响应函数的幅度和相位');