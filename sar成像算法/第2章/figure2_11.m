clear all;clc;close all;
A = figure();
N = 250;
len = {};
pslr = [];
wide_3dB = [];
for i = 0:1:6
    idx=num2str(i);
    figure(A);hold on;
    z = kaiser(N,i);
    plot(kaiser(N,i),'LineWidth',1);xlim([-10 260]);
    can = fft( z, 100*N ); 
    can = abs(can);
    can = 20.*log10( can./max(can) );
    can = fftshift( can );
    len = [len,['β=',idx]];
    pslr = [pslr,getpslr(can)];
    
    can1 = fft( z, 100*N );
    can1 = abs((can1).^2);
    can1 = 20.*log10( can1./max(can1) );
    can1 = fftshift( can1 );
    
    [t,t_r_location,t_l_location] = wide(can1);
    wide_3dB = [wide_3dB,t];


end

figure(A);hold on;
% plot(kaiser(N,2.5),'LineWidth',2,'color','red');xlim([-10 260]);
% len = [len,['β=2.5']];legend(len');
xlabel('kaiser窗的β'); title('图2.11 不同β值的kaiser窗形状');
% 画3dB宽度展宽比
q =  wide_3dB./wide_3dB(1)-1;
figure(2); subplot(1,2,1);
plot( 0:6,100*q);
xlabel('kaiser窗的β'); ylabel('展宽比'); title('3dB宽度展宽比')
grid on

% 画峰值旁瓣比
figure(2); subplot(1,2,2);
plot( 0:6, pslr );
xlabel('kaiser窗的β'); ylabel('PSLR,dB'); title('峰值旁瓣比(PSLR)')
grid on
sgtitle('图2.12 不同kaiser窗的展宽和峰值旁瓣比');

function [t,t_r_location,t_l_location] = wide(a)
[m,max_location] = max(a);
[r,t_r_location] = min(abs(a(max_location:end)+3));
[l,t_l_location] = min(abs(a(1:max_location)+3));
t_r_location = t_r_location + max_location - 1;
t = t_r_location - t_l_location;
end


function[k] = getpslr(a)
p = findpeaks(a);
p = sort(p,'descend');
k = p(2);
end
