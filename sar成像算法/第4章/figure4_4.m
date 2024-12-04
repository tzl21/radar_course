close all;clear all;clc
% 参数设置
Re = 6368e3; % 本地地球半径
h  = 800e3;  % 平台距地面高度
dr = 13.6  ; % 斜距分辨率

% 设置要扫描的所有地距
G_all = linspace( 255, 425, 101 )*1e3;
Ng    = size(G_all);

% 初始化结果
R0_all  = zeros( Ng ); % 斜距
thi_all = zeros( Ng ); % 入射角
dg_all  = zeros( Ng ); % 地距分辨率

%% 迭代
for cg=1:1:Ng(1,2)

    % 提取地距
    G = G_all(cg);

    % 计算地心角
    betae = G/Re;

    % 计算斜距R0
    temp = ( Re + h )^2 + Re^2 - 2*( Re + h )*Re*cos( betae );
    R0   = sqrt( temp );

    % 计算星下点离线角
    temp = Re/(R0/sin(betae));
    thn  = asin( temp );

    % 计算入射角
    thi = betae + thn;

    % 计算地距分辨率
    dg = dr/sin(thi);

    % 结果记录
    R0_all (cg)= R0;
    thi_all(cg)= thi;
    dg_all(cg) = dg;
end
%% 画图
figure,set(gcf,'Color','w');
subplot(3,1,1),plot(G_all/1e3, R0_all/1e3);
ylabel('斜距(km)');
subplot(3,1,2),plot(G_all/1e3, thi_all*180./pi);
ylabel('入射角(°)');
subplot(3,1,3),plot(G_all/1e3, dg_all);
xlabel('地距(km)'),ylabel('地距分辨率(m)');
sgtitle('图4.4 RADARSAT-1 W1波束地距分辨率的变化')