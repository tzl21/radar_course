clear;close all;clc;
u = [1,3,-1,5,2,6,4,-2];%信号
v = [1,2,3];%滤波器
w = conv(u,v);%卷积

figure;
subplot(3,1,1);stem(u);   %信号图
axis([-1 11 -7 26]);title('信号')
subplot(3,1,2);stem(v);   %滤波器图
axis([-1 11 -7 26]);title('滤波器')
subplot(3,1,3);stem(w);   %卷积图
axis([-1 11 -7 26]);title('卷积')
sgtitle('图2.1 线性卷积运算中的信号');