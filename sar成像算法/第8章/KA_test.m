%% 2023.12.2 wK仿真(Stolt插值)
clc;clear;close all
%% 参数设置
% 距离向参数

Tr = 40e-6;                      % 发射脉冲时宽
Kr = 0.5e+12;                    % 距离向调频率
% alpha_os_r = 1.2;               % 距离过采样率
% Nrg = 3072;                     % 距离采样点数
% % 距离向参数
Fr = 24e+6;
% Bw = abs(Kr)*Tr;                % 距离信号带宽
% Fr = alpha_os_r*Bw;             % 距离向采样率
% Nr = round(Fr*Tr);              % 距离采样点数
% 方位向参数
c = 299792458;                  % 电磁传播速度
Vr = 7100;                       % 等效雷达速度
Vs = 1.06*Vr;                        % 卫星平台速度
Vg = 0.94*Vr;                        % 波束扫描速度
f0 = 5.3e+9;                    % 雷达工作频率
Fa = 1700;
Ta = 0.64;
% alpha_os_a = 1.25;              % 方位过采样率
% Naz = 512;                     % 方位采样点数
Naz = round(2*Ta*Fa);
Nrg = round(2*Tr*Fr);
theta_sqc = +4*pi/180;
theta_r_c = asin(sin(theta_sqc)*Vr/Vg);         % 波束斜视角
R_ref = 1e+6;
R_eta_c = (R_ref + 20e+3)/cos(theta_r_c);                % 景中心斜距
% 方位向参数
lambda = c/f0;                  % 雷达工作波长
t_eta_c = -R_eta_c*sin(theta_r_c)/Vr;
                                % 波束中心穿越时刻
f_eta_c = 2*Vr*sin(theta_r_c)/lambda;
                                % 多普勒中心频率
      
                                % 目标照射时间
R0 = R_eta_c*cos(theta_r_c);    % 最短斜距
Ka = 2*Vr^2*cos(theta_r_c)^2/lambda/R0;              
La = 10;  
theta_bw = 0.886*lambda/La;     % 方位向3dB波束宽度
theta_syn = Vs/Vg*theta_bw;     % 合成角宽度
Ls = R_eta_c*theta_syn;         % 合成孔径长度
Delta_f_dop = 2*Vs*cos(theta_r_c)/c*f0*theta_bw; 
 
%  参数计算
rho_r = c/(2*Fr);               % 距离向分辨率
rho_a = La/2;                   % 距离向分辨率
Trg = Nrg/Fr;                   % 发射脉冲时宽
Taz = Naz/Fa;                   % 目标照射时间
d_t_tau = 1/Fr;                 % 距离采样时间间隔
d_t_eta = 1/Fa;                 % 方位采样时间间隔
d_f_tau = Fr/Nrg;               % 距离采样频率间隔    
d_f_eta = Fa/Naz;               % 方位采样频率间隔
%% 目标设置
% NPosition = [R0-1500*cos(theta_r_c),-1500*sin(theta_r_c);
%              R0-750*cos(theta_r_c),-750*sin(theta_r_c);
%              R0,0;
%              R0+750*cos(theta_r_c),+750*sin(theta_r_c);
%              R0+1500*cos(theta_r_c),+1500*sin(theta_r_c);]; 
NPosition = [R0,0];
% 设置数组
%  得到目标点的波束中心穿越时刻
[Ntarget,~] = size(NPosition);
Tar_t_eta_c = zeros(1,Ntarget);
for i = 1 :Ntarget
    DeltaX = NPosition(i,2) - NPosition(i,1)*tan(theta_r_c);
    Tar_t_eta_c(i) = DeltaX/Vs;
end
%  得到目标点的绝对零多普勒时刻
Tar_t_eta_0 = zeros(1,Ntarget);
for i = 1 : Ntarget
    Tar_t_eta_0(i) = NPosition(i,2)/Vr;
end
%% 变量设置
%  时间变量 以景中心的零多普勒时刻作为方位向零点
t_tau = (-Trg/2:d_t_tau:Trg/2-d_t_tau) + 2*R_eta_c/c;   % 距离时间变量
t_eta = (-Taz/2:d_t_eta:Taz/2-d_t_eta) + t_eta_c;       % 方位时间变量
%  长度变量
r_tau = (t_tau*c/2)*cos(theta_r_c);                     % 距离长度变量
%  频率变量 
f_tau = (-Fr/2:d_f_tau:Fr/2-d_f_tau);           % 距离频率变量 
f_eta = ((-Fa/2:d_f_eta:Fa/2-d_f_eta)+f_eta_c); % 方位频率变量
%% 坐标设置     
%  以距离时间为X轴，方位时间为Y轴
[t_tauX,t_etaY] = meshgrid(t_tau,t_eta);                % 设置距离时域-方位时域二维网络坐标
%  以距离长度为X轴，方位频率为Y轴                                                                                                            
[r_tauX,f_etaY] = meshgrid(r_tau,f_eta);                % 设置距离时域-方位频域二维网络坐标
%  以距离频率为X轴，方位频率为Y轴                                                                                                            
[f_tau_X,f_eta_Y] = meshgrid(f_tau,f_eta);              % 设置频率时域-方位频域二维网络坐标
%% 信号设置：原始回波信号                                                                                                  
st_tt = zeros(Naz,Nrg);
for i = 1 : Ntarget
    %  计算目标点的瞬时斜距
    R_eta = sqrt( NPosition(i,1)^2 +Vr^2*(t_etaY-Tar_t_eta_0(i)).^2 );   
    %  距离向包络
    wr = (abs(t_tauX-2*R_eta/c) <= Tr/2);                               
    %  方位向包络
    wa = sinc(0.886*atan(Vr*(t_etaY-t_eta_c)/R0)/theta_bw).^2;
    %  接收信号叠加
    st_tt_tar = wr.*wa.*exp(-1j*4*pi*f0*R_eta/c).*exp(+1j*pi*Kr*(t_tauX-2*R_eta/c).^2);                                                        
    st_tt = st_tt + st_tt_tar;  
end

% 绘图：原始数据
figure(1)                
subplot(221),imagesc( real(st_tt))             
xlabel('距离向(采样点)'),ylabel('方位向(采样点)'),title('(a)实部');
subplot(222),imagesc( imag(st_tt))
xlabel('距离向(采样点)'),ylabel('方位向(采样点)'),title('(b)虚部');
subplot(223),imagesc(  abs(st_tt)) 
xlabel('距离向(采样点)'),ylabel('方位向(采样点)'),title('(c)幅度');
subplot(224),imagesc(angle(st_tt))
xlabel('距离向(采样点)'),ylabel('方位向(采样点)'),title('(c)相位');


%% 信号设置：二维FFT
S0 = st_tt.*exp(-1i*2*pi*f_eta_c*t_etaY); %对齐到方位向多普勒中心；
S2df = fftshift(fft2(fftshift(S0)));
% S2df = (fft2((S0)));
% 绘图：二维fft
figure(2)                
subplot(111),imagesc( abs(S2df))             
xlabel('距离频率(采样点)'),ylabel('方位频率(采样点)'),title('二维FFT');

%% 信号设置:参考函数相乘（一致压缩）
theta_ref = 4*pi*R_ref / c * sqrt((f0+f_tau_X).^2-c^2*f_eta_Y.^2/(4*Vr^2)) + pi*f_tau_X.^2/Kr;
Href = exp(1j * theta_ref);
% 线性相位补偿
S_RFM = S2df .* Href.*exp(-1j*2*pi*(f_tau_X*(t_tau(Nrg/2+1)))).* exp(-1j*2*pi*(f_etaY*(t_eta(Naz/2+1))));
% 绘图：一致压缩
figure(3)                
subplot(111),imagesc( abs(S_RFM))             
xlabel('距离频率(采样点)'),ylabel('方位频率(采样点)'),title('一致压缩');
% 绘图：二维IFFT
figure(4)                
subplot(211),imagesc( abs(ifftshift(ifft2(ifftshift(S_RFM)))));             
xlabel('距离向(采样点)'),ylabel('方位向(采样点)'),title('二维IFFT（无stolt插值）');
subplot(212),imagesc( abs((ifft(ifftshift(S_RFM),Nrg,2))));             
xlabel('距离向(采样点)'),ylabel('方位向(采样点)'),title('RD（无stolt插值）');

%% 二维频域进行stolt插值
% 频率轴映射
f_new_c = sqrt((f0+0).^2-(c*(f_etaY)/Vr/2).^2) - f0;                        % 映射后距离向中心频率
fnew_m = f_tau_X + round((f_new_c - f_tau_X)/Fr)*Fr;                        % 将频率[-Fr/2, Fr/2]映射回其实际（卷绕前）对应的频率 
fold_m = sqrt((f0+fnew_m).^2+(c*(f_etaY)/Vr/2).^2)-f0;                      % 由新的频率，计算旧的频率


% 插值
Sstolt = zeros(Naz,Nrg);                                                    % 新建矩阵，存放插值后的数据
for ii = 1:Naz
    for jj = 1:Nrg
        Delta = (fold_m(ii,jj)-f_tau(jj))/(Fr/Nrg);                         % 计算相对频率差，并量化
        IntNum = ceil(Delta);                                               % 频率差向上取整
        kk = jj+IntNum;                                                     % 在原矩阵中的索引
        if(5<=kk && kk<=Nrg-3)                                              % 边界限定
            DecNum = IntNum-Delta;                                          % 频率差的小数部分
            SincVal = sinc(DecNum-4:1:DecNum+3)';                           % 生成插值核
            Sstolt(ii,jj) = S_RFM(ii,kk-4:kk+3)*SincVal;                    % 插值
        end
    end
end

%绘图：stolt插值
figure(5)                
subplot(111),imagesc( abs(Sstolt))             
xlabel('距离频率(采样点)'),ylabel('方位频率(采样点)'),title('stolt插值');

%% 信号设置：二维IFFT
s_image = ifftshift(ifft2(ifftshift(Sstolt))); 
% 绘图：二维IFFT
figure(6)                
subplot(111),imagesc( abs(s_image))             
xlabel('距离向(采样点)'),ylabel('方位向(采样点)'),title('二维IFFT（有stolt插值）');



%% 信号设置:点目标B的分析
for iiiii=1:1:Ntarget
    len = 32;%点附近len*len框
    cut = -len/2:1:len/2-1;

    ppp = round(Naz/2 +  NPosition(iiiii,2)/Vr*Fa);
    qqq = round(Nrg/2 +  2*(NPosition(iiiii,1)-R0)/c*Fr);
    start_tt = s_image(ppp-len:ppp+len,qqq-len:qqq+len);
    [p,q]=find(abs(start_tt)==max(max(abs(start_tt))));
    start_tt = s_image(ppp-(len-p)-1 + cut,qqq-(len-q)-1 + cut);
    
    Start_ff = fft2(start_tt);
 
    freq=16;%升采样倍数
    S_ac_test_2=fft2(start_tt);
    S_ac_test_buling_1 = zeros(len,freq*len);
    S_ac_test_buling   = zeros(freq*len,freq*len);
    for pp=1:1:len         
            [mi,I] = min(abs(S_ac_test_2(pp,:)));
            S_ac_test_buling_1(pp,1:I) = S_ac_test_2(pp,1:I);
            S_ac_test_buling_1(pp,freq*len-(len-I)+1:freq*len) = S_ac_test_2(pp,I+1:len);
    end
    for qq=1:1:freq*len        
            [mi,I]  = min(abs(S_ac_test_buling_1(:,qq)));
            S_ac_test_buling(1:I,qq) = S_ac_test_buling_1(1:I,qq);
            S_ac_test_buling(freq*len-(len-I)+1:freq*len,qq) = S_ac_test_buling_1(I+1:len,qq);
            for j=I+1:1:len
                 S_ac_test_buling(freq*len-(len-I)+j-I,qq) = S_ac_test_buling_1(j,qq);
             end
    end
    Dataq = S_ac_test_buling;
    dataq = ifft2(Dataq);
    zuihou =abs(dataq)/max(max(abs(dataq)));
    zuihou2 = 20*log10(zuihou);
    
    figure(7+iiiii);contour(abs(dataq),16);
    xlabel('距离时间(采样点)'),ylabel('方位时间(采样点)'),title('升采样后的点目标',iiiii)
  
    % 将两个方向的格点拉均匀
    [~,p_test]= max(abs(zuihou));
    [v_test,~] = max(abs(zuihou));
    [~,q_test] = max(v_test);
    i = p_test(q_test); %二维矩阵最大值，所在第几行
    j = q_test;      %二维矩阵最大值，所在第几列
    I=zuihou;
    I = I/max(max(abs(I)));
    [m,n] = size(I);
    % 变换后，两个方向上的格点距离均为juli_gedian
    fangwei_gedian = d_t_eta*Vs;
    juli_gedian = d_t_tau*c/2;
    m = round(m*fangwei_gedian/juli_gedian);  
    gedian = juli_gedian;
    I_new = imresize(I,[m n],'bicubic');
    I_new_log = 20*log10(abs(I_new)/max(max(abs(I_new))));

    if m<n
        if mod(m,2)==1
            m=m-1;
        end
      I_new2 = I_new(:,round(n/2)-m/2+1:1:round(n/2)+m/2);
    else
    if mod(n,2)==1
           n=n-1;
    end
     I_new2 = I_new(round(m/2)-n/2+1:1:round(m/2)+n/2,:);
    end
    I_new2_log = 20*log10(abs(I_new2)/max(max(abs(I_new2))));

    I_new3_real = imrotate(real(I_new2),theta_r_c/pi*180,'bicubic','crop');%,'crop'
    I_new3_imag = imrotate(imag(I_new2),theta_r_c/pi*180,'bicubic','crop');%,'crop'
    I_new3 = I_new3_real + 1j*I_new3_imag;
    I_new3_log = 20*log10(abs(I_new3)/max(max(abs(I_new3))));

    [m,n] = find(I_new3_log ==max(max(I_new3_log)));
    figure(7);
    subplot(2,Ntarget,iiiii+0*Ntarget),plot(I_new3_log(m,:)); title('点目标距离向剖面图');
    axis([-inf inf -40 0]);grid on;
    subplot(2,Ntarget,iiiii+1*Ntarget),plot(I_new3_log(:,n)); title('点目标方位向剖面图');
    axis([-inf inf -40 0]);grid on;

    [irw_a,locleft_3_a,locright_3_a] = get_irw(I_new3(:,n));
    [pslr_a] = get_pslr(I_new3_log(:,n));
    [islr_a] = get_islr(I_new3(:,n),2.25);

    [irw_r,locleft_3_r,locright_3_r] = get_irw(I_new3(m,:));
    [pslr_r] = get_pslr(I_new3_log(m,:));
    [islr_r] = get_islr(I_new3(m,:),2.25);

    disp("目标"+iiiii+"  距离分辨率："+irw_r/freq*c/2/Fr+"  方位分辨率："+irw_a/freq/Fa*Vr+"  距离PSLR:"+pslr_r+"  方位PSLR："+pslr_a+"  距离ISLR:"+islr_r+"  方位ISLR："+islr_a);
end

function [irw,locleft,locright] = get_irw(Af)
    % 找到Af的最大位置
    [~,locmax] = max(Af);
    % 找到locmax左边最接近-3dB的位置
    [~,locleft] = min(abs(Af(1:locmax-1)/max(abs(Af(1:locmax-1)))-1/sqrt(2)));
    % 找到locmax右边最接近-3dB的位置
    [~,locright] = min(abs(Af(locmax+1:end)/max(abs(Af(locmax+1:end)))-1/sqrt(2)));
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
function [islr] = get_islr(Af,alpha)%主瓣宽度为irw的alpha倍
    % 找到Af的最大位置
    [~,locmax] = max(Af);
    % 找到locmax左边最接近-3dB的位置
    [~,locleft] = min(abs(Af(1:locmax-1)/max(abs(Af(1:locmax-1)))-0.707));
    % 找到locmax右边最接近-3dB的位置
    [~,locright] = min(abs(Af(locmax+1:end)/max(abs(Af(locmax+1:end)))-0.707));
    locright = locright + locmax;
    % 得到3dB波束宽度
    irw = locright-locleft;
    locleft=locleft-round(alpha*irw/2);
    locright=locright+round(alpha*irw/2);
    % 计算总功率
    P_total = sum(Af(1:end).^2);
    % 计算主瓣功率
    P_main = sum(Af(locleft:locright).^2);
    % 一维积分旁瓣比
    islr = 10*log10((P_total-P_main)./P_main);
end
