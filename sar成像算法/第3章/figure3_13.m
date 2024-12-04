close all;clear all;clc;

T = 10e-6;          % 脉冲持续时间
B = 15e6;           % 脉冲带宽
K = B/T;            % 调频率
ratio = 5;          % 过采样率
Fs = ratio*B;       % 采样频率
dt = 1/Fs;          % 采样间隔
Nr = ceil(T/dt);    % 采样点数
t0 = ((0:Nr-1)-Nr/2)/Nr*T;          % 基本时间轴

st0 = exp(1i*pi*K*(t0-T/5).^2);     % 基本信号
space1 = zeros(1,round(Nr/5));      % 生成空信号
space2 = zeros(1,Nr);               % 生成空信号
st = [space1,st0,space2,st0,space2,st0,space1];     % 实际信号

N = length(st);             % 实际信号长度
n = 0:N-1;                  % 样本轴
f = ((0:N-1)-N/2)/N*Fs;     % 基本频率轴

Sf = fftshift(fft(st));                         % 实际信号的傅里叶变换
Hf1 = fftshift(fft(conj(fliplr(st0)),N));       % 方式1的匹配滤波器：时间反褶后取复共轭，计算N点补零DFT
Hf2 = fftshift(conj(fft(st0,N)));               % 方式2的匹配滤波器：补零后计算DFT，对结果取复共轭
Hf3 = exp(1i*pi*f.^2/K);                        % 方式3的匹配滤波器：直接在频域生成匹配滤波器

out1 = ifft(ifftshift(Sf.*Hf1));
out2 = ifft(ifftshift(Sf.*Hf2));
out3 = ifft(ifftshift(Sf.*Hf3));

figure,set(gcf,'Color','w');
subplot(4,1,1),plot(n,real(st));axis tight
title('(a)输入阵列信号的实部');ylabel('幅度');
subplot(4,1,2),plot(n,abs(out1));axis tight
title('(b)方式1的匹配滤波输出');ylabel('幅度');
subplot(4,1,3),plot(n,abs(out2));axis tight
title('(c)方式2的匹配滤波输出');ylabel('幅度');
subplot(4,1,4),plot(n,abs(out3));axis tight
title('(d)方式3的匹配滤波输出');xlabel('时间(采样点)');ylabel('幅度');
sgtitle('图3.13 通过压缩目标的位置来说明基带信号的弃置区和TA值');