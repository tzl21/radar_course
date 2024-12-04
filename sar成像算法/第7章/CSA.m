%% 2023.10.28 CSA仿真
clc
clear
close all
%% 参数设置
% 距离向参数
R_eta_c = 20e+3;                % 景中心斜距
Tr = 2.5e-6;                    % 发射脉冲时宽
Kr = 20e+12;                    % 距离向调频率
alpha_os_r = 1.6;               % 距离过采样率
Nrg = 2048;                      % 距离线采样点数
% 距离向参数
Bw = abs(Kr)*Tr;                % 距离信号带宽
Fr = alpha_os_r*Bw;             % 距离向采样率
Nr = round(Fr*Tr);              % 距离采样点数
% 方位向参数
c = 299792458;                  % 电磁传播速度
Vr = 150;                      % 等效雷达速度
Vs = Vr;                        % 卫星平台速度
Vg = Vr;                        % 波束扫描速度
f0 = 5.3e+9;                    % 雷达工作频率
Delta_f_dop = 120;             % 多普勒带宽
alpha_os_a = 1.2;              % 方位过采样率
Naz = 800;                      % 方位采样点数
theta_r_c = +8*pi/180;        % 波束斜视角
% 方位向参数
lambda = c/f0;                  % 雷达工作波长
t_eta_c = -R_eta_c*sin(theta_r_c)/Vr;
                                % 波束中心穿越时刻
f_eta_c = 2*Vr*sin(theta_r_c)/lambda;
                                % 多普勒中心频率
La = 0.886*2*Vs*cos(theta_r_c)/Delta_f_dop;               
                                % 实际天线长度
Fa = alpha_os_a*Delta_f_dop;    % 方位向采样率
Ta = 0.886*lambda*R_eta_c/(La*Vg*cos(theta_r_c));
                                % 目标照射时间
R0 = R_eta_c*cos(theta_r_c);    % 最短斜距
Ka = 2*Vr^2*cos(theta_r_c)^2/lambda/R0;              
                                % 方位向调频率
theta_bw = 0.886*lambda/La;     % 方位向3dB波束宽度
theta_syn = Vs/Vg*theta_bw;     % 合成角宽度
Ls = R_eta_c*theta_syn;         % 合成孔径长度
%  参数计算
rho_r = c/(2*Fr);               % 距离向分辨率
rho_a = La/2;                   % 方位向分辨率
Trg = Nrg/Fr;                   % 发射脉冲时宽
Taz = Naz/Fa;                   % 目标照射时间
d_t_tau = 1/Fr;                 % 距离采样时间间隔
d_t_eta = 1/Fa;                 % 方位采样时间间隔
d_f_tau = Fr/Nrg;               % 距离采样频率间隔    
d_f_eta = Fa/Naz;               % 方位采样频率间隔
%% 目标设置
NPosition = [R0-1e3,0;
             R0,0;
             R0+1e3,0
             R0+1e3,100];    % 设置点目标的（斜距，方位）
%  得到目标点的波束中心穿越时刻
Ntarget = length(NPosition);
Tar_t_eta_c = zeros(1,Ntarget);
for i = 1 : Ntarget
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
f_tau = fftshift(-Fr/2:d_f_tau:Fr/2-d_f_tau);           % 距离频率变量 
f_eta = fftshift((-Fa/2:d_f_eta:Fa/2-d_f_eta)+f_eta_c); % 方位频率变量

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
    R_eta = sqrt( NPosition(i,1)^2 +...
                  Vr^2*(t_etaY-Tar_t_eta_0(i)).^2 );   
    %  距离向包络
    wr = (abs(t_tauX-2*R_eta/c) <= Tr/2);                               
    %  方位向包络
    wa = (abs(t_etaY-Tar_t_eta_c(i))<= Ta/2);
    %  接收信号叠加
    st_tt_tar = wr.*wa.*exp(-1j*4*pi*f0*R_eta/c).*exp(+1j*pi*Kr*(t_tauX-2*R_eta/c).^2);                                                        
    st_tt = st_tt + st_tt_tar;  
end

% 绘图：原始数据
figure;imagesc( real(st_tt))             
xlabel('距离向(采样点)'),ylabel('方位向(采样点)'),title('原始信号实部');


%% 信号设置：变换到距离多普勒域，进行“补余RCMC”  
st_tt = st_tt.*exp(-1j*2*pi*f_eta_c*t_etaY);                    % 方位频率搬移到零多普勒中心（把f_eta_c弄没）
S_RD = fft(st_tt,Naz,1);                                        % 进行方位向傅里叶变换，得到距离多普勒域频谱
D_fn_Vr = sqrt(1-lambda^2.*(f_eta.').^2./(4*Vr^2));             % 徙动因子，列向量
D_fn_Vr_mtx = D_fn_Vr*ones(1,Nrg);                              % 形成矩阵
D_fn_ref_Vr = sqrt(1-lambda^2*f_eta_c^2/(4*Vr^2));              % 参考频率f_eta_c处的徙动因子，是常数
K_src = 2*Vr^2*f0^3.*D_fn_Vr.^3./(c*R0*(f_eta.').^2);           % 列向量，使用R0处的值 
K_src_mtx = K_src*ones(1,Nrg);                                  % 形成矩阵
Km = Kr./(1-Kr./K_src_mtx);                                     % 矩阵，这是变换到距离多普勒域的距离调频率
% 下面生成 变标方程 s_sc
s_sc = exp(1j*pi.*Km.*(D_fn_ref_Vr./D_fn_Vr_mtx-1).*(t_tauX-2*R0./(c.*D_fn_Vr_mtx)).^2);
% 下面将距离多普勒域的信号与变标方程相乘，实现“补余RCMC”
S_RD_1 = S_RD.*s_sc;       

% 绘图：方位向fft+补余RCMC
figure;imagesc( abs(S_RD))             
xlabel('距离向(采样点)'),ylabel('方位频率(采样点)'),title('距离多普勒域');
figure;imagesc( abs(S_RD_1))
xlabel('距离向(采样点)'),ylabel('方位频率(采样点)'),title('补余RCMC');

%% 信号设置：变换到二维频域，进行“距离压缩，SRC，一致RCMC”    
S_2df_1 = fft(S_RD_1,Nrg,2);                        % 进行距离向FFT，变换到二维频域。距离零频在两端
% 完成距离压缩，SRC，一致RCMC这三者相位补偿的滤波器H1：
H1 = exp(1j*pi.*D_fn_Vr_mtx./(D_fn_ref_Vr.*Km).*f_tau_X.^2)...
    .*exp(1j*4*pi/c.*(1./D_fn_Vr_mtx-1/D_fn_ref_Vr).*R0.*f_tau_X);
S_2df_2 = S_2df_1.*H1;    	                        % 在二维频域，相位相乘，实现距离压缩，SRC，一致RCMC
S_RD_2 = ifft(S_2df_2,Nrg,2);                       % 进行距离IFFT，回到距离多普勒域，完成所有距离处理。

% 绘图：方位向IFFT+距离压缩+SRC+一致RCMC
figure;imagesc(abs(S_RD_2))
xlabel('距离向(采样点)'),ylabel('方位频率(采样点)'),title('完成距离压缩，SRC，一致RCMC');

%% 信号设置：距离多普勒域，完成"方位压缩”和“附加相位校正”
R0_RCMC = (c/2).*t_tau*cos(theta_r_c);                         % 随距离线变化的R0，记为R0_RCMC
% 生成方位向匹配滤波器
Haz = exp(1j*4*pi.*(D_fn_Vr*R0_RCMC).*f0./c);                  % 方位MF
Haf_offset = exp(-1j*2*pi*f_etaY.*t_eta_c);                    % 使景中心成像到图像的方位向中心!!!方位向存在中心频率f_eta_c
% 附加相位校正项
H2 = exp(-1j*4*pi.*Km./(c^2).*(1-D_fn_Vr_mtx./D_fn_ref_Vr)...
    .*((1./D_fn_Vr)*R0_RCMC-R0./D_fn_Vr_mtx).^2); 	           % 附加相位校正项
H3 = exp(j*4*pi*R0/lambda);                                    % 恢复目标的全部相位
% 下面进行相位相乘，在距离多普勒域，同时完成方位MF和附加相位校正
S_RD_3 = S_RD_2.*Haz.*H2.*Haf_offset;                          % 距离多普勒域，相位相乘
% 最后通过IFFT回到图像域，完成方位处理
s_image = ifft(S_RD_3,Naz,1).*H3; 	                           % 完成成像过程，H3为恢复全部的相位

% 绘图：方位压缩+相位校正
figure;imagesc(abs(S_RD_3))             
xlabel('距离向(采样点)'),ylabel('方位频率(采样点)'),title('方位压缩和附加相位校正');
figure;imagesc(abs(s_image))
xlabel('距离向(采样点)'),ylabel('方位向(采样点)'),title('成像结果');


%% 信号设置:点目标B的分析
for iiiii=1:1:Ntarget
    len = 32;%点附近len*len框
    cut = -len/2:1:len/2-1;
    n = round(Naz/2 +  NPosition(iiiii,2)/Vr*Fa);
    m = round(Nrg/2 +  2*(NPosition(iiiii,1)-R0)/c*Fr);
    start_tt = s_image(n-len:n+len,m-len:m+len);
    [p,q]=find(abs(start_tt)==max(max(abs(start_tt))));
    start_tt = s_image(n-(len-p)-1 + cut,m-(len-q)-1 + cut);
    
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
    
    figure(7);
    subplot(3,Ntarget,iiiii),imagesc(abs(dataq))
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
    subplot(3,Ntarget,iiiii+1*Ntarget),plot(I_new3_log(m,:)); title('点目标距离向剖面图');
    axis([-inf inf -40 0]);grid on;
    subplot(3,Ntarget,iiiii+2*Ntarget),plot(I_new3_log(:,n)); title('点目标方位向剖面图');
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
