clear all;clc;close all;
N=4;
sn=2^5;
sf=2^8;
t =linspace( -N/2, N/2, sn );
s =sinc(t);
% 得到fs
dt = (t(length(t))-t(1))/(length(t)-1);
fs = 1/dt;
% 得到频谱; *dt转化为FT的结果
af =fftshift( fft( s, sf) )* dt;
% 得到频点
fre = (0:sf-1)*fs/sf - fs/2;

figure
subplot( 2, 2, 1); 
plot(fre, abs(af) );xlim([-1,1]);ylim([0,1.2]); xlabel('频率');ylabel('幅度');title('(a)4点的sinc插值核');

N=8;
sn=2^5;
sf=2^8;
t =linspace( -N/2, N/2, sn );
s =sinc(t);
% 得到fs
dt = (t(length(t))-t(1))/(length(t)-1);
fs = 1/dt;
% 得到频谱; *dt转化为FT的结果
af =fftshift( fft( s, sf) ) * dt;
% 得到频点
fre = (0:sf-1)*fs/sf - fs/2;

subplot( 2, 2, 2); 
plot(fre, abs(af) );xlim([-1,1]); ylim([0,1.2]);xlabel('频率');ylabel('幅度');title('(b)8点的sinc插值核');

N=16;
sn=2^5;
sf=2^8;
t =linspace( -N/2, N/2, sn );
s =sinc(t);
% 得到fs
dt = (t(length(t))-t(1))/(length(t)-1);
fs = 1/dt;
% 得到频谱; *dt转化为FT的结果
af =fftshift( fft( s, sf) ) * dt;
% 得到频点
fre = (0:sf-1)*fs/sf - fs/2;

subplot( 2, 2, 3); 
plot(fre, abs(af) );xlim([-1,1]); ylim([0,1.2]);xlabel('频率');ylabel('幅度');title('(c)16点的sinc插值核');

N=8;
sn=2^5;
sf=2^8;
t =linspace( -N/2, N/2, sn );
s =sinc(t);
%加窗
beta=2.5;
s = s .* kaiser( sn, beta )';
% 得到fs
dt = (t(length(t))-t(1))/(length(t)-1);
fs = 1/dt;
% 得到频谱; *dt转化为FT的结果
af =fftshift( fft( s, sf) ) * dt;
% 得到频点
fre = (0:sf-1)*fs/sf - fs/2;

subplot( 2, 2, 4); 
plot(fre, abs(af) );xlim([-1,1]);ylim([0,1.2]); xlabel('频率');ylabel('幅度');title('(d)加ksiser窗的8点的sinc插值核');

sgtitle('图2.16 sinc函数插值核的频谱');
