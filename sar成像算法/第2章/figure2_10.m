clear all;clc;close all;
fs = 400;
F_continue = 0:0.1:2*fs;
F_apparent_complex = F_continue-(round(F_continue./fs)).*fs;
F_apparent_real = abs(F_apparent_complex);
subplot(2,1,1);plot(F_continue,F_apparent_complex);
title('(a)复信号');xlabel('可观测频率');ylabel('可观测频率');
set(gca,'XTick',0:200:800);set(gca,'XTicklabel',{'0','0.5fs','fs','1.5fs','2fs'})
subplot(2,1,2);plot(F_continue,F_apparent_real);
title('(b)实信号');xlabel('可观测频率');ylabel('可观测频率');
set(gca,'XTick',0:200:800);set(gca,'XTicklabel',{'0','0.5fs','fs','1.5fs','2fs'})
sgtitle('图2.10 由混叠引起的可观测频率');