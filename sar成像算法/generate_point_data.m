function [  s0, f_etac, delta_r, delta_a, center_R0, center_etac  ] = generate_point_data( parameters )
%GENERATE_POINT_DATA generate SAR target raw data
%   parameters is a struct (including kinds of simulation parameters as its
%   properties.
%   Tips:
%       Please give a complete parameters struct.
%       If parameters.NUM_TARGETS equals to zero, this function will
%       display a scene, guiding you to set proper targets.
%       
%   例如：
%         NUM_TARGETS = 3;    % 仿真的目标数为3
%         rs = [0, 0, 20];    % 各目标距离向距离
%         as = [-20, 0, -10]; % 目标相对方位向距离
%         parameters = struct(...
%             'center_Rc', center_Rc,...          % 景中心斜距
%             'theta_rc_deg', theta_rc_deg,...    % 斜视角
%             'Nrg', Nrg,...                      % 距离向采样点数
%             'Naz', Naz,...                      % 方位向采样点数
%             'Vr', Vr,...                        % 载机速度
%             'f0', f0,...                        % 载波频率
%             'Tr', Tr,...                        % 发射脉冲宽度
%             'Kr', Kr,...                        % 发射脉冲调频率
%             'BW_dop', BW_dop,...                % 多普勒带宽
%             'alpha_os_r', alpha_os_r,...        % 距离向过采样率
%             'alpha_os_a', alpha_os_a,...        % 方位向过采样率
%             'NUM_TARGETS', NUM_TARGETS,...      % 点目标数量
%             'rs', rs,...                        % 点目标距离向坐标（m）
%             'as', as...                         % 点目标方位向坐标（m）
%         );
%         [ s0, f_etac, delta_r, delta_a, center_R0, center_Rc ] = generate_point_data(parameters);
%
%   返回值:
%   s0 仿真回波信号矩阵
%   f_etac 信号多普勒中心频率
%   delta_r 数据距离向采样间距（m）
%   delta_a 数据方位采样间距（m）
%   center_R0 景中心点最近斜距
%   center_Rc 景中心波束中心穿越时刻斜距

% parse the parameters
center_Rc = parameters.center_Rc;   % 景中心点波束穿越时刻斜距
theta_rc = parameters.theta_rc_deg * pi / 180;
Nrg = parameters.Nrg;
Naz = parameters.Naz;
Vr = parameters.Vr;
f0 = parameters.f0;
Tr = parameters.Tr;
Kr = parameters.Kr;
BW_dop = parameters.BW_dop;
alpha_os_r = parameters.alpha_os_r;
alpha_os_a = parameters.alpha_os_a;

% targets info
NUM_TARGETS = parameters.NUM_TARGETS;
rs = parameters.rs;
as = parameters.as;

% default parameters
c = 299792458;    % light speed

% derived parameters
Vs = Vr;
Vg = Vr;
lambda = c / f0;
center_R0 = center_Rc * cos(theta_rc);          % 景中心点最短斜距
center_etac = - center_Rc * sin(theta_rc) / Vr; % 景中心点对应的相对波束中心穿越时刻;
La = 0.886 * 2 * Vs * cos(theta_rc) / BW_dop;   % 天线孔径长度
beta_bw = 0.886 * lambda / La;                  % we suppose
Fr = Tr * Kr * alpha_os_r;  % 距离采样率
Fa = BW_dop * alpha_os_a;   % 脉冲重复频率
delta_r = c/2/Fr;           % 距离向采样间距（SAR信号空间）
delta_a = Vr / Fa;          % 方位向采样间距
% not used by generating the signal, and they will be returned directly for
% the user
f_etac = 2 * Vr * sin(theta_rc) / lambda; % 多普勒中心频率

% 观测时间轴区间确定
tau0 = 2 * center_Rc / c;
eta0 = center_etac;
tau = ((-Nrg / 2) : (Nrg / 2 - 1)) / Fr + tau0;
eta = ((-Naz / 2 : Naz / 2 - 1)) / Fa + eta0;

%##########################################################################
% 在此探究可以设置的点目标坐标取值范围
% SAR原始信号空间的距离向和方位向时间宽度
tau_w = tau(end) - tau(1);
eta_h = eta(end) - eta(1);
% 允许的斜距坐标范围
% 用linspace有问题，仍需优化，参考extra/better_generate_time_freq_axis.m
r_w = linspace(-(tau_w/2-Tr/2), (tau_w/2-Tr/2), Nrg) * c / 2 * cos(theta_rc); % -Tp/2是因为考虑了脉冲宽度
% 斜视角使得不同斜距处目标点允许的方位位置不同(1)
a_h_top = eta_h / 2 * Vg + r_w * tan(theta_rc);
a_h_bottom = -eta_h / 2 * Vg + r_w * tan(theta_rc);
% % 拓宽eta轴，以便最终能成像容纳所有点（不在操作（2）之后进行是因为eta轴还要足够宽确保原始信号能被完全接收）
% % 因为user已经定好了Naz，我们就不要擅自改变了
% Naz = round((a_h_top(end)-a_h_bottom(1)) / interval_a);
% Naz = ceil(Naz/2)*2;
% % disp(['Naz2:', num2str(Naz)]);
% eta = ((-Naz / 2 : Naz / 2 - 1)) / Fa + eta0;
% 斜距的不同使得不同斜距处目标时域（r-a）能量轨迹方位向长度不同(2)
a_h_top = a_h_top - (r_w + center_R0) / cos(theta_rc) * beta_bw / 2;
a_h_bottom = a_h_bottom + (r_w + center_R0) / cos(theta_rc) * beta_bw / 2;

% 显示能放置的目标点的范围
figure;
plot(r_w, a_h_top, 'r-');
hold on;
plot(r_w, a_h_bottom, 'r-');
line([r_w(1); r_w(1)], [a_h_bottom(1); a_h_top(1)], 'color', 'r');
line([r_w(end); r_w(end)], [a_h_bottom(end); a_h_top(end)], 'color', 'r');
xlabel('相对景中心斜距向r（m）');ylabel('方位向（m）');
set(gca, 'YDir', 'reverse');
grid on;

if NUM_TARGETS == 0
    % just want to know the coordiante range of target we can set
    % suptitle('斜距-方位平面下目标点可设置场景坐标范围示意图');
    s0 = [];
    hold off;
    return;
else
    % suptitle('斜距-方位平面下目标点设置场景示意图');
    plot(rs, as, '*');
    hold off;
end

%##########################################################################
% 计算各目标点的最近斜距
R0s = center_R0 + rs;
% 计算各目标点的绝对零多普勒时刻
eta_0s = as / Vr;
% 计算各目标点的绝对多普勒中心时刻
eta_cs = eta_0s + center_etac - rs * tan(theta_rc) / Vr;

%% 3. 构造雷达原始数据
% 时间坐标网格
[tau_mtx, eta_mtx] = meshgrid(tau, eta);

% 计算目标点瞬时斜距（随方位时间变化）
% 注：etaY-eta_0s(i)为目标点相对各自零多普勒时刻的方位向时间eta
R_eta = zeros(NUM_TARGETS, Naz, Nrg);
for i = 1:NUM_TARGETS
    R_eta(i, :, :) = sqrt((R0s(i)^2 + Vr^2 * (eta_mtx - eta_0s(i)).^2 ));
end

A0 = 1;
s0 = zeros(Naz, Nrg);
for i = 1:NUM_TARGETS
    % 包络
    w_r = (abs(tau_mtx - 2 * reshape(R_eta(i, :, :), Naz, Nrg) / c) <= Tr / 2);
    w_a = sinc(0.886 / beta_bw * atan(Vg * (eta_mtx - eta_cs(i)) / R0s(i))).^2;
    % 相位
    theta1 = -4 * pi * f0 * reshape(R_eta(i, :, :), Naz, Nrg) / c;
    theta2 = pi * Kr * (tau_mtx - 2 * reshape(R_eta(i, :, :), Naz, Nrg) / c).^2;
    % 信号多点累加
    s0 = s0 + A0 * w_r .* w_a .* exp(1j*theta1) .* exp(1j*theta2);
end
end

