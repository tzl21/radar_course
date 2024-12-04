clear all;clc;close all;
t = -100:0.01:100-0.01;
fs = 200;
A = zeros(1,20000);
for i = 9500:1:10499
    A(i) = 1;
end
subplot(2,2,1)
plot(t,A)
title('时间t');
xlim([-40,40])
N=length(t);
Af = fftshift(fft(fftshift(A)));
f=-fs/2:fs/N:(fs/2-fs/N);
subplot(2,2,2);plot(f,Af);title('频率f');xlim([-4,4])

dt = -100:0.01:100-0.01;
fs = 200;
N=length(dt);
sincx=sinc(dt);
subplot(2,2,3)
plot(dt,sincx);
title('时间t');
xlim([-4,4])
sincf = fftshift(fft(fftshift(sincx)));
df=-fs/2:fs/N:(fs/2-fs/N);
subplot(2,2,4);plot(df,sincf);title('频率f');xlim([-2,2]);
sgtitle('图2.3 矩形和sinc函数的傅里叶变换');