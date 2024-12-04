close all;clear all;clc
%信号持续时间T=7.24us，信号带宽B=5.8MHz，将过采样率设为5是为了更清晰地观测信号波形。

T = 7.24e-6;                % 信号持续时间
B = 5.8e6;                  % 信号带宽
K = B/T;                    % 调频率
ratio = 5;                  % 过采样率
Fs = ratio*B;               % 采样频率
dt = 1/Fs;                  % 采样间隔
N = ceil(T/dt);             % 采样点数
t=-T/2:T/N:T/2-T/N;         % 时间轴

st = exp(1i*pi*K*t.^2);     % 生成信号

figure,set(gcf,'Color','w');
subplot(2,2,1),plot(t*1e6,real(st));
title('(a)信号的实部'),xlabel('相对于t_{0}时间(\mus)'),ylabel('幅度');
subplot(2,2,2),plot(t*1e6,pi*K*t.^2);
title('(c)信号的相位'),xlabel('相对于t_{0}时间(\mus)'),ylabel('弧度');
subplot(2,2,3),plot(t*1e6,imag(st));
title('(b)信号的虚部'),xlabel('相对于t_{0}时间(\mus)'),ylabel('幅度');
subplot(2,2,4),plot(t*1e6,K*t*1e-6);
title('(d)信号频率'),xlabel('相对于t_{0}时间(\mus)'),ylabel('MHz');
sgtitle('图3.1 线性调频信号的相位和频率');