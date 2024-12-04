clc
clear
close all

% 设置参数TBP和T
TBP = 100;
T   = 10e-6;

% 求解B和K
B = TBP/T;
K = B/T;

% 生成时间序列t
t = linspace( -T, T, 2001 );

% 计算第一部分结果p1 第二部分结果p2 完整结果a
% 第一部分对T归一化了
p1 = ( T - abs(t) ) ./ T;
p2 = sinc( K.*t.* ( T-abs(t) ) );
a  = p1.*p2;

% % 图片设置
% h = figure; set( h, 'position', [100,20,600,650]);
% sub_row = 3; sub_col = 1; sub_count = 0;
% % 画图 慢变部分
% sub_count = sub_count + 1;  subplot( sub_row, sub_col, sub_count );
% plot( t*1e6, p1 ); ylabel('第一部分'); title('匹配滤波器输出慢变部分')
% % 画图 快变部分
% sub_count = sub_count + 1;  subplot( sub_row, sub_col, sub_count );
% plot( t*1e6, p2 ); ylabel('第二部分'); title('匹配滤波器输出快变部分')
% % 画图 慢变部分
% sub_count = sub_count + 1;  subplot( sub_row, sub_col, sub_count );
% plot( t*1e6, a ); ylabel('两部分'); title('匹配滤波器输出'); xlabel('时间（us）')

figure
plot( t*1e6, a ); ylabel('幅度'); xlabel('时间（us）');
title('图3.5 匹配滤波器输出的3dB分辨率的测量');
