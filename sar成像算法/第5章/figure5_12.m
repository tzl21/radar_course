clc
clear
close all
% 距离向参数
R_eta_c = 850e+3;                       % 景中心斜距
% 方位向参数
Vr = 7100;                              % 等效雷达速度
Ta = 0.64;                              % 目标照射时间
Ka = 2095;                              % 方位向调频率
theta_r_c = -0.3*pi/180;                % 斜视角
% 参数计算
Delta_f_dop = abs(Ka*Ta);               % 方位信号带宽
t_eta_c = -R_eta_c*sin(theta_r_c)/Vr;   % 波束中心穿越时刻
% 参数设置
alpha_a_s = 1;                          % 方位过采样率
Fa = alpha_a_s*Delta_f_dop;             % 方位采样频率PRF
Na = 2*ceil(Fa*Ta/2);                   % 方位采样点数
dt = Ta/Na;                             % 采样时间间隔
df = Fa/Na;                             % 采样频率间隔 
% 变量设置
t_eta =  (-Ta/2:dt:Ta/2-dt) + t_eta_c;  % 方位时间变量
% 信号表达
R_eta = R_eta_c -...
        Vr*sin(theta_r_c)*(t_eta-t_eta_c) +...
        (1/2)*(Vr^2*cos(theta_r_c)^2/R_eta_c)*(t_eta-t_eta_c).^2;   % 瞬时斜距展开式
RCM_1 = -Vr*sin(theta_r_c)*(t_eta-t_eta_c+Ta/2);                    % 距离徙动线性分量
RCM_2 = (1/2)*(Vr^2*cos(theta_r_c)^2/R_eta_c)*(t_eta-t_eta_c).^2;   % 距离徙动二次分量
RCM_all = RCM_1 + RCM_2;                                            % 距离徙动总量
% 绘图  
figure(1),set(gcf,'Color','w');
plot(RCM_1,t_eta,'k--');
hold on;
plot(RCM_all,t_eta,'r');
xlabel('距离徙动(m)'),ylabel('方位时间(s)');axis([-5 35 0 1.2]);
sgtitle('图5.12 距离单元徙动的线性和二次分量')