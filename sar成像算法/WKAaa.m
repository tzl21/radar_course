close all;clear all;
%% 1. 仿真参数 (参考 p142, table 6.1)
center_Rc = 20e3;  % 景中心斜距
Vr = 150;       % 等效雷达速度
Tr = 2.5e-6;    % 发射脉冲时宽
Kr = 40e12;     % 距离调频率
f0 = 5.3e9;     % 雷达工作频率
BW_dop = 80;    % 多普勒带宽
Fr = 120e6;  % 距离采样率
Fa = 100;   % 方位采样率
Naz = 256;  % 方位向采样点数（距离线条数）
Nrg = 3072;  % 距离向采样点数（距离线采样点数）
theta_rc_deg = 5; % 低斜视角5度
c = 299792458;    % 光速

% derived params
lambda = c / f0;
theta_rc = theta_rc_deg * pi / 180;
Vs = Vr;
Vg = Vr;
Np = Tr * Fr;   % 脉冲序列长度（采样点数）
alpha_os_r = Fr / (Kr*Tr);
alpha_os_a = Fa / BW_dop;

%% 2. 生成原始雷达数据
NUM_TARGETS = 7;    % 仿真的目标数
gap = 50 * 8;
% gap = 25;
rs = [-3, -2, -1, 0, 0, 2, 3]*gap;
as = [-3, -2, -1, 0, 1, 2, 3]*gap*tan(theta_rc);
parameters = struct(...
    'center_Rc', center_Rc,...          % 景中心斜距
    'theta_rc_deg', theta_rc_deg,...    % 斜视角
    'Nrg', Nrg,...                      % 距离向采样点数
    'Naz', Naz,...                      % 方位向采样点数
    'Vr', Vr,...                        % 载机速度
    'f0', f0,...                        % 载波频率
    'Tr', Tr,...                        % 发射脉冲宽度
    'Kr', Kr,...                        % 发射脉冲调频率
    'BW_dop', BW_dop,...                % 多普勒带宽
    'alpha_os_r', alpha_os_r,...        % 距离向过采样率
    'alpha_os_a', alpha_os_a,...        % 方位向过采样率
    'NUM_TARGETS', NUM_TARGETS,...      % 点目标数量
    'rs', rs,...                        % 点目标距离向坐标（m）
    'as', as...                         % 点目标方位向坐标（m）
);

[ s0, f_etac, delta_r, delta_a, center_R0, center_Rc ] = generate_point_data(parameters);

%% 3. 成像处理
Rref = center_R0;
s2 = wKA( s0, lambda, Kr, Vr, Fr, Fa, Rref, f_etac, Tr );


%% show
Vg = Vr;
[Naz, Nrg] = size(s2);
x = ((-Nrg / 2) : (Nrg / 2 - 1)) / Fr * c / 2;
y = ((-Naz / 2 : Naz / 2 - 1)) / Fa * Vg;
figure;
imagesc(x, y, abs(s2));
xlabel('距离向（m）');ylabel('方位向（m）');
title('完整处理后的目标');set(gca, 'YDir', 'normal');

%% 4. 点目标分析
% 计算每个点出现位置的索引值
ns = round(rs/delta_r) + (Nrg/2 + 1);
ms = round(as/delta_a) + (Naz/2 + 1);
len = 16;

figure; % 放大显示每个点目标
for i = 1:NUM_TARGETS
    target = s2(ms(i)-len/2:ms(i)+len/2-1, ns(i)-len/2:ns(i)+len/2-1);
    subplot(ceil(NUM_TARGETS/3), 3, i);
    imagesc(abs(target));
    xlabel('距离向（采样点）');ylabel('方位向（采样点）');ylabel('方位向（采样点）');title(['目标', num2str(i)]);
end

% 升采样分析点目标A（即第一个点目标）
p = 1;
target = s2(ms(p)-len/2:ms(p)+len/2-1, ns(p)-len/2:ns(p)+len/2-1);
[image_upsample,signal_r,quality_r,signal_a,quality_a] = f_point_analyse(target,delta_r,delta_a);

BW_r= abs(Kr*Tr);
La = 0.886 * 2 * Vs * cos(theta_rc) / BW_dop;   % 天线孔径长度
IRW_r_theory = c/2/BW_r*0.886*1.18;
IRW_a_theory = La/2*Vg/Vs*1.185;
disp(['距离向理论分辨率:',num2str(IRW_r_theory),'m']);
disp(['方位向理论分辨率:',num2str(IRW_a_theory),'m']);

