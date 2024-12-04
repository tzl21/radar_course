clc
clear
close all
% 距离向参数
R_eta_c = 450;                              % 景中心斜距
% 方位向参数
c  = 3e+8;                                  % 电磁传播速度
f0 = 5.3e+6;                                % 雷达工作频率
Vr = 7100;                                  % 等效雷达速度
Ta = 0.64;                                  % 目标照射时间
Ka = 2095;                                  % 方位向调频率
theta_r_c = 0*pi/180;                       % 斜视角
% 参数计算
lambda = c/f0;                              % 雷达工作波长
R0 = R_eta_c*cos(theta_r_c);                % 最短斜距
Delta_f_dop = abs(Ka*Ta);                   % 方位信号带宽
% 参数设置
alpha_a_s = 1;                              % 方位过采样率
Fa = alpha_a_s*Delta_f_dop;                 % 方位采样频率PRF
Na = 2*ceil(Fa*Ta/2);                       % 方位采样点数
dt = Ta/Na;                                 % 采样时间间隔
df = Fa/Na;                                 % 采样频率间隔 
% 变量设置
t_eta_1 =  (-Ta/16:dt:Ta/16-dt);            % 方位时间变量
t_eta_2 =  (-Ta/16:dt:Ta/16-dt)-0.2;        % 方位时间变量
t_eta_3 =  (-Ta/16:dt:Ta/16-dt)-0.5;        % 方位时间变量
f_eta =  (-Fa/16:df:Fa/16-df);              % 方位频率变量
% 信号表达
R_eta_1 = R0 + Vr^2*t_eta_1.^2/(2*R0);      % 目标轨迹
R_eta_2 = R0 + Vr^2*(t_eta_2+0.2).^2/(2*R0);% 目标轨迹
R_eta_3 = R0 + Vr^2*(t_eta_3+0.5).^2/(2*R0);% 目标轨迹
% 信号表达                                          
R_rd_1 = R0 + lambda^2*R0/(8*Vr^2)*f_eta.^2;    % 目标轨迹
R_rd_2 = R0 + lambda^2*R0/(8*Vr^2)*f_eta.^2;    % 目标轨迹
R_rd_3 = R0 + lambda^2*R0/(8*Vr^2)*f_eta.^2;    % 目标轨迹
% 绘图  
figure(1),set(gcf,'Color','w');
subplot(1,3,1);
plot(R_eta_1,t_eta_1),hold on
plot(R_eta_2,t_eta_2),hold on
plot(R_eta_3,t_eta_3)
title('RA域');xlabel('距离向'),ylabel('方位向');
subplot(1,3,2)
plot(f_eta,t_eta_1),hold on
plot(f_eta,t_eta_2),hold on
plot(f_eta,t_eta_3)
title('DA域');xlabel('方位频率'),ylabel('方位向');
subplot(1,3,3);
plot(R_eta_1,f_eta),hold on
plot(R_eta_2,f_eta),hold on
plot(R_eta_3,f_eta)
title('RD域');xlabel('距离向'),ylabel('方位频率');
sgtitle('图5.14 在距离多普勒域中，多个目标轨迹重合为一个轨迹')