close all;clear all;clc
Ta = 64;                      % 脉冲持续时间
Ka = -1.56e-2;                   % 方位向调频率
% 参数计算
Delta_f_dop = abs(Ka*Ta);      % 方位信号带宽
% 参数设置
ratio =1;                      % 方位过采样率
Fa = ratio*Delta_f_dop;        % 方位采样频率PRF
N = 2*ceil(Fa*Ta/2);           % 方位采样点数
dt = Ta/N;                     % 采样时间间隔
df = Fa/N;                     % 采样频率间隔
% 变量设置
t1 = -Ta/2:dt:Ta/2;            % 时间变量
t2 = -3*Ta/4:dt:Ta/4;          % 时间变量
f1 = -Fa/2:df:Fa/2;            % 频率变量
f2 = -3*Fa/4:df:Fa/4;          % 频率变量
n=0:1:N;
% 信号表达
st1 = exp(1j*pi*Ka*t1.^2);     % Chirp信号复数表达式
st2 = exp(1j*pi*Ka*t2.^2);     % Chirp信号复数表达式
omega_a = sinc(1.5*t1./Ta).^2; % 双程天线波束方向图
st1  = omega_a.*st1;            % 天线方向图调制的方位Chirp信号
st2  = omega_a.*st2;            % 天线方向图调制的方位Chirp信号
omega_a_nor = abs(omega_a)/max(abs(omega_a));       % 归一化
omega_a_log = 20*log10(omega_a_nor);                % 对数化
Sf1 = fft(st1);
Sf2 = fft(st2);
% 绘图
figure(1),set(gcf,'Color','w');
subplot(2,2,1),plot(n,real(st1));
title('(a)信号实部(斜视角为零)');ylabel('幅度');axis([0,N-1,-1.2,1.2]);
subplot(2,2,2),plot(n,abs(Sf1));
title('(b)信号频谱(斜视角为零)');ylabel('幅度');axis([0,N-1,0,10]);
subplot(2,2,3),plot(n,real(st2));
title('(c)信号实部(斜视角非零)');xlabel('时间(采样点)');ylabel('幅度');axis([0,N-1,-1.2,1.2]);
subplot(2,2,4),plot(n,abs(Sf2));
title('(d)信号频谱(斜视角非零)');xlabel('频率(采样点)');ylabel('幅度');axis([0,N-1,0,10]);
sgtitle('图5.5 斜视角为零和非零时的多普勒中心')