clear all;clc;close all;
%% 参数
T = 1;
dt = 1/400;
t = 0:dt:1;
fs = 1/dt;
f0 = 15;
ft = 100;
BW = 75;%带宽
k = BW/T;%调频率
N  = length(t);

%% 幅度修正
p0 = t+3;
p_max = max(p0);
p = p0./p_max;

%% 建立原信号
s1 = cos(2*pi*(f0.*t+0.5*k*t.^2)).*p;
s2 = exp(1i*2*pi*((f0-52.5).*t+0.5*k*t.^2)).*p;
s3 = cos(2*pi*(ft.*t+0.5*k*t.^2)).*p;
s4 = exp(1i*2*pi*(ft.*t+0.5*k*t.^2)).*p;

%% 傅里叶变换
F1_1 = fftshift(fft(s1,N));
F1 = F1_1./max(abs(F1_1));
F2_2 = fftshift(fft(s2,N));
F2 = F2_2./max(abs(F2_2));
F3_3 = fftshift(fft(s3,N));
F3 = F3_3./max(abs(F3_3));
F4_4 = fftshift(fft(s4,N));
F4 = F4_4./max(abs(F4_4));

%% 绘图
f = -fs/2:fs/N:fs/2-fs/N;
figure
subplot(2,2,1); plot( f,abs(F1));ylabel('幅度'); title('(a)实&基带');grid on
subplot(2,2,2); plot( f,abs(F2));ylabel('幅度'); title('(b)复&基带');grid on
subplot(2,2,3); plot( f,abs(F3));ylabel('幅度'); title('(c)实&非基带'); grid on
subplot(2,2,4); plot( f,abs(F4));ylabel('幅度'); title('(d)复&非基带'); grid on
sgtitle('图2.6 实信号和复信号，基带信号和非基带信号间的频谱差异');