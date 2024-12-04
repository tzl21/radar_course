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
t1 = -Ta/2:dt1:Ta/2-dt1;           % 时间变量
t2 = -Ta/2:dt2:Ta/2-dt2;           % 时间变量
f1 = -Fa1/2:df1:Fa1/2-df1;         % 频率变量
f2 = -Fa2/2:df2:Fa2/2-df2;         % 频率变量
% 信号表达
st1 = exp(1j*pi*Ka*t1.^2);     % Chirp信号复数表达式
st2 = exp(1j*pi*Ka*t2.^2);     % Chirp信号复数表达式
omega_a = sinc(6*t2./Ta).^2;    % 双程天线波束方向图
st  = omega_a.*st2;            % 天线方向图调制的方位Chirp信号
omega_a_nor = abs(omega_a)./max(abs(omega_a));       % 归一化
omega_a_log = 20*log10(omega_a_nor);                % 对数化

% 窗函数
window_1 = kaiser(N2,2.5);                         % 时域窗
Window_1 = fftshift(window_1);                      % 频域窗
% 信号变换-->方式一
ht_1 = conj(fliplr(st));                            % 将时间反褶后的复制脉冲取复共轭
ht_window_1 = window_1'.*ht_1;                       % 加窗
Hf_1 = fftshift(fft(ht_window_1,N2));               % 计算补零离散傅里叶变换
% 窗函数
window_2 = kaiser(N2,2.5);                         % 时域窗
Window_2 = fftshift(window_2);                      % 频域窗
% 信号变换-->方式二
ht_2 = st;                                          % 复制信号
ht_window_2 = window_2'.*ht_2;                       % 加窗
Hf_2 = fftshift(conj(fft(ht_window_2,N2)));         % 计算补零离散傅里叶变换
% 窗函数
window_3 = kaiser(N2,2.5)';                         % 时域窗
Window_3 = fftshift(window_3);                      % 频域窗
% 信号变换-->方式三
Hf_3 = exp(1j*pi*f2.^2/Ka);                         % 计算补零离散傅里叶变换
% 信号表达
Sf = fftshift(fft(st));
Sf_1 = Sf.*Hf_1;
st_1 = ifft(ifftshift(Sf_1));                       % 方式一匹配滤波结果
st_1_nor = abs(st_1)./max(abs(st_1));               % 归一化
st_1_log = 20*log10(st_1_nor);                      % 对数化
Sf_2 = Sf.*Hf_2;
st_2 = ifft(ifftshift(Sf_2));                       % 方式二匹配滤波结果
st_2_nor = abs(st_2)./max(abs(st_2));               % 归一化
st_2_log = 20*log10(st_2_nor);                      % 对数化
Sf_3 = Sf.*Hf_3;
st_3 = ifft(ifftshift(Sf_3));                       % 方式三匹配滤波结果
st_3_nor = abs(st_3)./max(abs(st_3));               % 归一化
st_3_log = 20*log10(st_3_nor);                      % 对数化

% 绘图
figure(1),set(gcf,'Color','w');
subplot(3,1,1),plot(t2,real(st2));
title('(a)方位Chirp信号实部，忽略天线波束方向图');ylabel('幅度');axis([-Ta/2-5,Ta/2+5,-1.2,1.2]);
subplot(3,1,2),plot(t2,omega_a_log);
title('(b)双程天线波束方向图');ylabel('幅度(dB)');axis([-Ta/2-5,Ta/2+5,-60,0]);
subplot(3,1,3),plot(t2,fftshift(st_1_log));
title('(c)压缩后的目标与虚影(鬼影)');xlabel('方位时间');ylabel('幅度(dB)');axis([-Ta/2-5,Ta/2+5,-50,0]);
sgtitle('图5.4 由方位chirp信号混叠造成的方位模糊')
