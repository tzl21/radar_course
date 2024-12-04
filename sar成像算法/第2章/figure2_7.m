clear all;clc;close all;
T = 10e-3;
f = 300;
t = 0:T/10000:T-T/10000;
S = sin(2*pi*f*t);

f1 = 800;
d1 = 1/f1;
t1 = 0:d1:T-d1;
f2 = 400;
d2 = 1/f2;
t2 = 0:d2:T-d2;
S1 = sin(2*pi*f*t1);

S2 = sin(2*pi*f*t2);
S2 = spline(t2,S2,t);%使用spline 基于非等间距的样本点对正弦曲线插值

%% 绘制图像
plot(t, S, 'b', t1, S1, '*',t, S2,'--');xlabel('时间 (s)');ylabel('幅度');
legend('原始sin函数(300Hz)', '800Hz采样点','400Hz采样点');
sgtitle('图2.7 以两个不同的采样率对300hz对正弦波采样来说明混叠现象');
