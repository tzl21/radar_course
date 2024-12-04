clear
% 时间序列生成
t = -8:0.02:8;
% 生成未经加权插值核
kernel1 = sinc( t );
% 经过加权插值核
kernel2 = sinc( t ) .* kaiser( length(t), 2.5 )';
% 作图
figure
plot( t, kernel1, 'b', t, kernel2, 'r');xlabel('时间（采样点）');ylabel('幅度');legend('未经加权','加权');
title('图2.15 Kaiser窗加权后的sinc函数，β=2.5');