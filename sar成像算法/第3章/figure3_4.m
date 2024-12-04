close all;clear all;clc
%信号持续时间T = 10us，信号带宽B = 10MHz，改变过采样率（1.4, 1.2, 1.0, 0.8）的取值以探索过采样率对信号频谱的影响。

T = 10e-6;                      % 信号持续时间
B = 10e6;                       % 信号带宽
K = B/T;                        % 调频率
ratio = [1.4,1.2,1.0,0.8];
Num = length(ratio);
figure,set(gcf,'Color','w');
for ii = 1:Num
    Fs = ratio(ii)*B;           % 采样频率
    dt = 1/Fs;                  % 采样间隔
    N = ceil(T/dt);             % 采样点数

    t = ((0:N-1)-N/2)/N*T;      % 时间轴
    f = ((0:N-1)-N/2)/N*Fs;

    st0 = exp(1i*pi*K*t.^2);    % 生成信号
    st1 = [st0,zeros(1,N)];     % 补零后的信号，补1倍
    Sf = fft(fftshift(st1));

    n = (0:2*N-1)/2;

    subplot(Num,2,2*ii-1),plot(t*1e6,real(st0));axis tight
    ylabel(['\alpha_{os}=',num2str(ratio(ii))]);
    if(ii == 1)
        title('信号实部');
    end
    if(ii == Num)
        xlabel('时间(\mus)');
    end
    subplot(Num,2,2*ii),plot(n,abs(Sf));axis tight
    if(ii ==1)
        title('频谱幅度');
    end
    if(ii == Num)
        xlabel('频率（单元）');
    end
end
sgtitle('图3.4 过采样率在频谱中引起的能量间隙');