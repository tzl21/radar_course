close all;clear all;clc
Ta = 128;                      % 脉冲持续时间
Ka = -0.085;                   % 方位向调频率
% 参数计算
Delta_f_dop = abs(Ka*Ta);      % 方位信号带宽
% 参数设置
ratio = [5,0.25];              % 方位过采样率
Fa1 = ratio(1)*Delta_f_dop;    % 方位采样频率PRF
N1 = 2*ceil(Fa1*Ta/2);         % 方位采样点数
dt1 = Ta/N1;                   % 采样时间间隔
df1 = Fa1/N1;                  % 采样频率间隔
Fa2 = ratio(2)*Delta_f_dop;    % 方位采样频率PRF
N2 = 2*ceil(Fa2*Ta/2);         % 方位采样点数
dt2 = Ta/N2;                   % 采样时间间隔
df2 = Fa2/N2;                  % 采样频率间隔
% 变量设置
t1 = -Ta/2:dt1:Ta/2;           % 时间变量
t2 = -Ta/2:dt2:Ta/2;           % 时间变量
f1 = -Fa1/2:df1:Fa1/2;         % 频率变量
f2 = -Fa2/2:df2:Fa2/2;         % 频率变量
% 信号表达
st1 = exp(1j*pi*Ka*t1.^2);     % Chirp信号复数表达式
st2 = exp(1j*pi*Ka*t2.^2);     % Chirp信号复数表达式
%参数计算
F1 = Ka*t1./Fa2;
F2 = (Ka*t2+floor((Fa2/2-Ka*t2)/Fa2)*Fa2)./Fa2;

% 绘图
figure,set(gcf,'Color','w');
subplot(3,1,1),plot(t1,real(st1));
title('(a)混叠前方位chirp信号实部');ylabel('幅度');axis([-Ta/2-5,Ta/2+5,-1.2,1.2]);
subplot(3,1,2),plot(t2,real(st2));
title('(b)混叠后方位chirp信号实部');ylabel('幅度');axis([-Ta/2-5 Ta/2+5,-1.2 1.2]);
subplot(3,1,3),plot(t1,F1,'b',t2,F2,'r');
title('(c)信号瞬时频率');xlabel('方位时间');ylabel('频率(PRF)');axis([-Ta/2-5 Ta/2+5,-2 2]);
sgtitle('图5.3 离散脉冲对方位信号的采样造成的方位混叠')

