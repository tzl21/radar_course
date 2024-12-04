clc
clear
close all

% 定义N 和 M
N = 2^3;
M = 16;

% 计算fs
Ts = 1;
fs = 1/Ts;

% 计算delta 记为d
d = 1/M;

% 生成数据n和m
n = 0:N-1;
m = (0:M)';
% 生成表格数据kerl
temp = N/2 - 1 + repmat( m*d, 1, N ) - repmat( n, M + 1, 1 );
kerl = sinc( temp ) ;
% 归一化
kerl = diag( 1./sqrt( sum(kerl.^2, 2) ) ) * kerl;

% 显示kerl 与table进行对比
disp('插值核为:\n');
disp( kerl(end:-1:1,:));


% 对kerl进行测试
% 生成抽样数据
t  = 0:N-1;
gn = sin( 2*pi*0.25*t );
% 生成画图数据
lt  = 0:0.01:N-1;
lgt = sin( 2*pi*0.25*lt ); 
% 插值
interp_data = kerl * gn';
interp_t    = N/2-1+m*d;

% 画图
figure
plot(t,gn,'ro'); hold on;
plot(lt,lgt);
plot( interp_t, interp_data, 'go')
grid on
legend('抽样数据','连续信号数据','插值数据')


