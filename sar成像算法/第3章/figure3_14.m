clc
clear
close all

% 参数设置
TBP = 42;                % 时间带宽乘积
T = 7.2e-6;              % 脉冲持续时间
N_st = 2^11;                         
% 参数计算
B = TBP/T;               % 信号带宽
K = B/T;                 % 调频频率
alpha_os = 200;          % 过采样率，使用过高的过采样率是为了方便地实现升采样
F = alpha_os*B;          % 采样频率
N = 2*ceil(F*T/2);       % 采样点数
dt = T/N;                % 采样时间间隔
df = F/N;                % 采样频率间隔
% 参数设置
t_c = 0e-6;              % 时间偏移
% 参数计算
f_c = -K*t_c;            % 中心频点
% 变量设置
t1 = -T/2:dt:T/2-dt;     % 时间变量
f1 = -F/2:df:F/2-df;     % 频率变量
% 信号表达                             
st = exp(1j*pi*K*(t1-t_c).^2);          % 发射信号
% % 绘图
% H1 = figure;
% set(H1,'position',[100,100,600,300]);
% subplot(211),plot(t1,real(st),'k')
% subplot(212),plot(t1,imag(st),'k')
% suptitle('发射信号')
% 参数设置
t_0 = 0e-6;              % 回波时延
% 变量设置
t2 = -T/2+t_0:dt:T/2+t_0-dt;            % 时间变量
f2 = -F/2+f_c:df:F/2+f_c-df;            % 频率变量
% 信号表达                                                                 
srt = exp(1j*pi*K*(t2-t_c-t_0).^2);     % 回波信号
% 信号表达
Srf = fftshift(fft(srt));

%%
% 参数设置
QPE = linspace(0,0.8*pi,N);            % 二次相位误差
dk = QPE/(pi*(T/2)^2);                 % 调频率误差
% 参数设置
IRW1  = zeros(1,N);                    % 初始化冲激响应宽度
PSLR1 = zeros(1,N);                    % 初始化峰值旁瓣比
ISLR1 = zeros(1,N);                    % 初始化积分旁瓣比
IRW2  = zeros(1,N);                    % 初始化冲激响应宽度
PSLR2 = zeros(1,N);                    % 初始化峰值旁瓣比
ISLR2 = zeros(1,N);                    % 初始化积分旁瓣比
IRW3  = zeros(1,N);                    % 初始化冲激响应宽度
PSLR3 = zeros(1,N);                    % 初始化峰值旁瓣比
ISLR3 = zeros(1,N);                    % 初始化积分旁瓣比
% 循环计算
for i = 1:N
    % 变量设置
    B_dk = (K+dk(i))*T;
    F_dk = alpha_os*B_dk;
    df_dk = F_dk/N;
    f3 = -F_dk/2+f_c:df_dk:F_dk/2+f_c-df_dk;           % 频率变量
    % 信号表达                                                                 
    st_dk = exp(1j*pi*(K+dk(i))*t1.^2);                  
    Sf_dk = fftshift(fft(st_dk));
    % 信号变换-->频域方式一
    window_1 = kaiser(N,2.5)';                         % 时域窗
    Window_1 = fftshift(window_1);                     % 频域窗
    ht_dk_1 = conj(fliplr(st_dk));                     % 将时间反褶后的复制脉冲取复共轭
    ht_window_dk_1 = window_1.*ht_dk_1;                % 加窗
    Hf_dk_1 = fftshift(fft(ht_window_dk_1,N));         % 计算补零离散傅里叶变换
    % 信号变换-->频域方式二
    window_2 = kaiser(N,2.5)';                         % 时域窗
    Window_2 = fftshift(window_2);                     % 频域窗
    ht_dk_2 = st_dk;                                   % 复制信号
    ht_window_dk_2 = window_2.*ht_dk_2;                % 加窗
    Hf_dk_2 = fftshift(conj(fft(ht_window_dk_2,N)));   % 计算补零离散傅里叶变换
    % 信号变换-->频域方式三
    window_3 = kaiser(N,2.5)';                         % 时域窗
    Window_3 = fftshift(window_3);                     % 频域窗
    Hf_dk_3 = Window_3.*exp(1j*pi*f3.^2/(K+dk(i)));    % 计算补零离散傅里叶变换
% 参数计算-->方式一                                         
    Soutf_dk_1 = Srf.*Hf_dk_1;
    soutt_dk_1 = ifft(ifftshift(Soutf_dk_1));          % 方式一匹配滤波结果 
    soutt_dk_1_nor = abs(soutt_dk_1)./max(abs(soutt_dk_1));               % 归一化
    soutt_dk_1_log = 20*log10(abs(soutt_dk_1)./max(abs(soutt_dk_1))+eps); % 对数化
    % 参数计算-->IRW
    [irw_dk_1,~,~] = get_irw(fftshift(soutt_dk_1_nor));
    IRW1(i) = irw_dk_1;
    % 参数计算-->PSLR
    [pslr_dk_1] = get_pslr(fftshift(soutt_dk_1_log));
    PSLR1(i) = pslr_dk_1;
    % 参数计算-->ISLR
    [islr_dk_1] = get_islr(fftshift(soutt_dk_1_nor),5);
    ISLR1(i) = islr_dk_1;
% 参数计算-->方式二                                        
    Soutf_dk_2 = Srf.*Hf_dk_2;
    soutt_dk_2 = ifft(ifftshift(Soutf_dk_2));          % 方式二匹配滤波结果 
    soutt_dk_2_nor = abs(soutt_dk_2)./max(abs(soutt_dk_2));               % 归一化
    soutt_dk_2_log = 20*log10(abs(soutt_dk_2)./max(abs(soutt_dk_2))+eps); % 对数化
    % 参数计算-->IRW
    [irw_dk_2,~,~] = get_irw(fftshift(soutt_dk_2_nor));
    IRW2(i) = irw_dk_2;
    % 参数计算-->PSLR
    [pslr_dk_2] = get_pslr(fftshift(soutt_dk_2_log));
    PSLR2(i) = pslr_dk_2;
    % 参数计算-->ISLR
    [islr_dk_2] = get_islr(fftshift(soutt_dk_2_nor),5);
    ISLR2(i) = islr_dk_2;
% 参数计算-->方式三                                         
    Soutf_dk_3 = Srf.*Hf_dk_3;
    soutt_dk_3 = ifft(ifftshift(Soutf_dk_3));          % 方式三匹配滤波结果 
    soutt_dk_3_nor = abs(soutt_dk_3)./max(abs(soutt_dk_3));               % 归一化
    soutt_dk_3_log = 20*log10(abs(soutt_dk_3)./max(abs(soutt_dk_3))+eps); % 对数化
    % 参数计算-->IRW
    [irw_dk_3,~,~] = get_irw(soutt_dk_3_nor);
    IRW3(i) = irw_dk_3;
    % 参数计算-->PSLR
    [pslr_dk_3] = get_pslr(soutt_dk_3_log);
    PSLR3(i) = pslr_dk_3;
    % 参数计算-->ISLR
    [islr_dk_3] = get_islr(soutt_dk_3_nor,5);
    ISLR3(i) = islr_dk_3;
end
% 绘图
figure;
subplot(131),plot(QPE/pi,(IRW1-IRW1(1))/IRW1(1)*100)
title('(a)IRW'),xlabel('|QPE|(\pi弧度)'),ylabel('展宽百分比')
subplot(132),plot(QPE/pi,PSLR1)
title('(b)PSLR'),xlabel('|QPE|(\pi弧度)'),ylabel('PSLR/dB')
subplot(133),plot(QPE/pi,ISLR1)
title('(c)ISLR'),xlabel('|QPE|(\pi弧度)'),ylabel('ISLR/dB')
% subplot(334),plot(QPE/pi,(IRW2-IRW2(1))/IRW2(1)*100,'k')
% title('(a)IRW'),xlabel('|QPE|(\pi弧度)'),ylabel('展宽百分比')
% subplot(335),plot(QPE/pi,PSLR2,'k')
% title('(b)PSLR'),xlabel('|QPE|(\pi弧度)'),ylabel('PSLR/dB')
% subplot(336),plot(QPE/pi,ISLR2,'k')
% title('(c)ISLR'),xlabel('|QPE|(\pi弧度)'),ylabel('ISLR/dB')
% subplot(337),plot(QPE/pi,(IRW3-IRW3(1))/IRW3(1)*100,'k')
% title('(a)IRW'),xlabel('|QPE|(\pi弧度)'),ylabel('展宽百分比')
% subplot(338),plot(QPE/pi,PSLR3,'k')
% title('(b)PSLR'),xlabel('|QPE|(\pi弧度)'),ylabel('PSLR/dB')
% subplot(339),plot(QPE/pi,ISLR3,'k')
% title('(c)ISLR'),xlabel('|QPE|(\pi弧度)'),ylabel('ISLR/dB')
sgtitle('图3.14 当\beta=2.5时的IRW、PSLR、ISLR与二次相位误差之间的关系')



%% 提取冲击响应宽度
function [irw,locleft,locright] = get_irw(Af)
    % 找到Af的最大位置
    [~,locmax] = max(Af);
    % 找到locmax左边最接近-3dB的位置
    [~,locleft] = min(abs(Af(1:locmax-1)/max(abs(Af(1:locmax-1)))-0.707));
    % 找到locmax右边最接近-3dB的位置
    [~,locright] = min(abs(Af(locmax+1:end)/max(abs(Af(locmax+1:end)))-0.707));
    locright = locright + locmax;
    % 得到3dB波束宽度
    irw = locright-locleft;
end
%% 提取峰值旁瓣比
function [pslr] = get_pslr(Af)
    % 找到所有的pesks
    peaks = findpeaks(Af);
    % 对peaks进行排序
    peaks = sort(peaks,'descend');
    % 得到第一旁瓣
    pslr = peaks(2);
end
%% 提取积分旁瓣比
function [islr] = get_islr(Af,Nr)
    % 找到Af的最大位置
    [~,locmax] = max(Af);
    % 找到locmax左边最接近-3dB的位置
    [~,locleft] = min(abs(Af(1:locmax-1)/max(abs(Af(1:locmax-1)))-0.707));
    % 找到locmax右边最接近-3dB的位置
    [~,locright] = min(abs(Af(locmax+1:end)/max(abs(Af(locmax+1:end)))-0.707));
    locright = locright + locmax;
    % 计算总功率
    P_total = sum(Af(locleft-Nr:locright+Nr).^2);
    % 计算主瓣功率
    P_main = sum(Af(locleft:locright).^2);
    % 一维积分旁瓣比
    islr = 10*log10((P_total-P_main)./P_main);
end