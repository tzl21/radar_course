clear all;clc;close all;
% 设置时间长度T; 带宽BW; 初始频率f0
T  = 1;
BW = 300;
k  = BW/T;% 斜率
h = figure();
%% 初始复信号
% 生成时间
dt_0 = 1/1200;
t_0  = 0:dt_0:T;
fs_0 = 1/dt_0;
% 初始信号长度
N_0  = length(t_0);
% 幅度修正
p0_0 = t_0+2;
p_max_0 = max(p0_0);
p_0 = p0_0./p_max_0;
% 频率
f_0 = (0:N_0-1)*fs_0/N_0 - fs_0/2;

%% 采样复信号
% 生成时间
dt = 1/400;
t  = 0:dt:T;
fs = 1/dt;
% 初始信号长度
N  = length(t);
% 幅度修正
p0 = t+2;
p_max = max(p0);
p = p0./p_max;
% 频率
f = (0:N-1)*fs/N - fs/2;

%% 生成信号
for i = 1:5
    f0 = 30+(i-1)*100;
    %初始信号
    s0 = exp( 1i*2*pi*((f0-180).*t_0+0.5*k*t_0.^2)).*p_0;  % 复数基带
    F0 = fftshift(fft(s0,N_0));
    F0 = F0./max(abs(F0));
    %采样信号
    s = exp( 1i*2*pi*((f0-180).*t+0.5*k*t.^2)).*p;  % 复数基带
    F = fftshift(fft(s,N));
    F = F./max(abs(F));
    %绘图
    figure(h);
    subplot(5,2,2*(i-1)+1);plot(f_0,abs(F0));ylabel('幅度');grid on;xlim([-210,610])
    if i == 1, title('连续信号的的频谱');end
    subplot(5,3,3*(i-1)+3 );plot(f,abs(F));ylabel('幅度');grid on;
    if i == 1, title('采样后的频谱');end
end
sgtitle('图2.9 采样引起的频谱平移（复信号）');