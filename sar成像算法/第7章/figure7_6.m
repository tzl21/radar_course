clc, clear, close all;
%% 参数设置
%  雷达参数
c = 3e8;
Vr = 150;
f0 = 5.3e9;
%  距离参数
R0 = 20e3;
Fr = 60e6;
Naz = 256;
%  方位参数
Fa = 100;
Nrg = 5120;
f_eta_c = 0;
%% 参数计算
%  雷达参数
lambda = freq2wavelen(f0);
%% 变量设置
R_mtx = [R0-Nrg/4/Fr*c; R0; R0+Nrg/4/Fr*c];
f_eta = -Fa/2 : Fa/Naz : Fa/2-Fa/Naz;
%% 信号生成
RCM_total = lambda^2/8/Vr^2*R_mtx.*(f_eta.^2 - f_eta_c^2);
RCM_bulk  = lambda^2/8/Vr^2*R0*ones(3,1).*(f_eta.^2 - f_eta_c^2);
RCM_diff  = lambda^2/8/Vr^2*(R_mtx - R0*ones(3,1)).*(f_eta.^2 - f_eta_c^2);
%% 绘制图形
figure(1);
subplot(311), plot(R_mtx(1,:)/1000 + RCM_total(1,:), f_eta, ...
              R_mtx(2,:)/1000 + RCM_total(2,:), f_eta, ...
              R_mtx(3,:)/1000 + RCM_total(3,:), f_eta),
axis([10 30, -55 +55]), xlabel('距离/km'), ylabel('多普勒频率'), title('(a)整体距离徙动');
subplot(312), plot(R_mtx(1,:)/1000 + RCM_bulk(1,:), f_eta, ...
              R_mtx(2,:)/1000 + RCM_bulk(2,:), f_eta, ...
              R_mtx(3,:)/1000 + RCM_bulk(3,:), f_eta), 
axis([10 30, -55 +55]), xlabel('距离/km'), ylabel('多普勒频率'), title('(b)一致距离徙动');
subplot(313), plot(R_mtx(1,:)/1000 + RCM_diff(1,:), f_eta, ...
              R_mtx(2,:)/1000 + RCM_diff(2,:), f_eta, ...
              R_mtx(3,:)/1000 + RCM_diff(3,:), f_eta), 
axis([10 30, -55 +55]), xlabel('距离/km'), ylabel('多普勒频率'), title('(c)补余距离徙动');
sgtitle('图7.6 整体RCM可以表述为不随距离变化的一致RCM与随距离变化的补余RCM之和')