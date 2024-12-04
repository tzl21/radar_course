clear all;clc;close all;
%% 原始信号
T  = 1;%设置时间长度T
BW = 150;% 带宽BW
k  = BW/T;%斜率
h=figure();%使所有图在一张图中生成
dt_0 = 1/1200; 
t_0  = 0:dt_0:T;% 生成原始时间
fs_0 = 1/dt_0;
N_0 = length(t_0);
% 信号幅度修正
p0_0 = t_0+2;
p0_max_0 = max(max(p0_0));
p_0 = p0_0./p0_max_0;  
%% 抽样
dt = 1/400;
t  = 0:dt:T;
fs = 1/dt;
% 信号幅度修正
p0 = t+2;
p0_max = max(max(p0));
p = p0./p0_max;
N  = length(t);
% 频率
f = -fs/2:fs/N:fs/2-fs/N;
%% 构造信号与傅里叶变换
for i = 1:9
    f0 = 25+(i-1)*50;%初始频率f0
    %初始信号
    s_0 = cos(2*pi*(f0.*t_0+0.5*k*t_0.^2)).*p_0;
    % 进行fft
    F_0 = fftshift(fft(s_0,N_0));
    F_0 = F_0./max(abs(F_0));
    % 频率
    f_0 = (0:N_0-1)*fs_0/N_0 - fs_0/2;
    %采样信号
    s = cos(2*pi*(f0.*t+0.5*k*t.^2)).*p;
    % 进行fft
    F = fftshift(fft(s)./ N );  
    F = F./max(abs(F));
    % 画图
    figure(h); hold on
    subplot(9,2,2*(i-1)+1 );plot( f_0, abs(F_0) ); ylabel('幅度');grid on;
    title( ['FC=',num2str(f0+BW/2)])
    subplot(9,3,3*(i-1)+3 );plot( f, abs(F) ); ylabel('幅度');grid on;
    if i == 1
        title('采样后的频率');
    end
end
sgtitle('图2.8 采样引起的频谱平移（实信号）');

