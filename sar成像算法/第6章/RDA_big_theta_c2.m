%%  距离多普勒算法（RDA)大斜视角+方式2的二次压缩（精确）
%方式2：通过二维频域中的相位相乘予以实现。该域不能同时对距离时间和距离频率进行调整，
%故需假定R_0不变。如果需要,可以使用距离恒定区对R_0进行更新。
clc
clear
close all
%% 参数设置
% 距离向参数
R_eta_c = 20e+3;                % 景中心斜距
Tr = 2.5e-6;                    % 发射脉冲时宽
Kr = 20e+12;                    % 距离向调频率
alpha_os_r = 1.2;               % 距离过采样率
Nrg = 320;                      % 距离线采样点数
% 距离向参数
Bw = abs(Kr)*Tr;                % 距离信号带宽
Fr = alpha_os_r*Bw;             % 距离向采样率
Nr = round(Fr*Tr);              % 距离采样点数(脉冲序列长度)
% 方位向参数
c = 3e+8;                       % 电磁传播速度
Vr = 150;                       % 等效雷达速度
Vs = Vr;                        % 卫星平台速度
Vg = Vr;                        % 波束扫描速度
f0 = 5.3e+9;                    % 雷达工作频率
Delta_f_dop = 100;               % 多普勒带宽
alpha_os_a = 1.25;              % 方位过采样率
Naz = 256;                      % 距离线数
theta_r_c = +21.9*pi/180;       % 波束斜视角
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
rho_a = La/2;                   % 距离向分辨率
Trg = Nrg/Fr;                   % 发射脉冲时宽
Taz = Naz/Fa;                   % 目标照射时间
d_t_tau = 1/Fr;                 % 距离采样时间间隔
d_t_eta = 1/Fa;                 % 方位采样时间间隔
d_f_tau = Fr/Nrg;               % 距离采样频率间隔    
d_f_eta = Fa/Naz;               % 方位采样频率间隔
%% 目标设置
%  设置目标点相对于景中心之间的距离
A_r = -50; A_a = -50;                                   % A点位置
B_r = 0; B_a = +0;                                   % B点位置
C_r = +50; C_a = +50;        % C点位置
%  得到目标点相对于景中心的位置坐标
A_x = R0 + A_r; A_Y = A_a;                              % A点坐标
B_x = R0 + B_r; B_Y = B_a;                              % B点坐标
C_x = R0 + C_r; C_Y = C_a;                              % C点坐标
NPosition = [A_x,A_Y;
             B_x,B_Y;
             C_x,C_Y;];                                 % 设置数组
fprintf( 'A点坐标为[%+3.3f，%+3.3f]km\n', NPosition(1,1)/1e3, NPosition(1,2)/1e3 );
fprintf( 'B点坐标为[%+3.3f，%+3.3f]km\n', NPosition(2,1)/1e3, NPosition(2,2)/1e3 );
fprintf( 'C点坐标为[%+3.3f，%+3.3f]km\n', NPosition(3,1)/1e3, NPosition(3,2)/1e3 );
%  得到目标点的波束中心穿越时刻
Ntarget = 3;
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
f_tau = f_tau - round((f_tau-0)/Fr)*Fr;                 % 距离频率变量(可观测频率)                          
f_eta = fftshift(-Fa/2:d_f_eta:Fa/2-d_f_eta);           % 方位频率变量
f_eta = f_eta - round((f_eta-f_eta_c)/Fa)*Fa;           % 方位频率变量(可观测频率)
%% 坐标设置  
%  以距离时间为X轴，方位时间为Y轴
[t_tauX,t_etaY] = meshgrid(t_tau,t_eta);                % 设置距离时域-方位时域二维网络坐标
%  以距离长度为X轴，方位频率为Y轴                                                                                                            
[r_tauX,f_etaY] = meshgrid(r_tau,f_eta);                % 设置距离时域-方位频域二维网络坐标
%  以距离频率为X轴，方位频率为Y轴                                                                                                            
[f_tau_X,f_eta_Y] = meshgrid(f_tau,f_eta);              % 设置频率时域-方位频域二维网络坐标
%% 信号设置：原始回波信号                                                        
tic
wait_title = waitbar(0,'开始生成雷达原始回波数据 ...');  
pause(1);
st_tt = zeros(Naz,Nrg);
for i = 1 : Ntarget
    %  计算目标点的瞬时斜距
    R_eta = sqrt( NPosition(i,1)^2 +...
                  Vr^2*(t_etaY-Tar_t_eta_0(i)).^2 ); 
    %  后向散射系数幅度
    A0 = [1,1,1,1]*exp(+1j*0);   
    %  距离向包络
    wr = (abs(t_tauX-2*R_eta/c) <= Tr/2);                               
    %  方位向包络
    wa = sinc(0.886*atan(Vg*(t_etaY-Tar_t_eta_c(i))/NPosition(i,1))/theta_bw).^2;      
    %  接收信号叠加
    st_tt_tar = A0(i)*wr.*wa.*exp(-1j*4*pi*f0*R_eta/c)...
                            .*exp(+1j*pi*Kr*(t_tauX-2*R_eta/c).^2);                                                          
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
xlabel('距离向(采样点)'),ylabel('方位向(采样点)'),title('(d)相位');
sgtitle('图6.16 大斜视角（21.9°）时，三点目标的仿真原始数据')

pause(1);
close(wait_title);
toc

%% 信号设置：距离压缩
%  一次距离压缩
window = kaiser(Nrg,2.5)';                              % 时域窗
Window = fftshift(window);                              % 频域窗
%  计算滤波器
Hrf = Window.*exp(+1j*pi*f_tau_X.^2/Kr);  
%  匹配滤波
Sf_ft = fft(st_tt,Nrg,2).*Hrf;

% 绘图：不进行IFFT的距离压缩后的数据
figure(2)
subplot(221),imagesc(abs(Sf_ft))
ylabel('方位向(采样点)'),title('(a)距离压缩')


%% 信号设置--》方位向傅里叶变换
Saf_tf = fft(Sf_ft,Naz,1);
%  二次距离压缩
Drd0 = sqrt(1-lambda^2.*f_eta.^2/(4*Vr^2));             % 距离多普勒域中的徙动因子%%%%%%%%精确的是f_eta为变量
Ksrc = 2*Vr^2*f0^3.*Drd0.^3/(c*R0.*f_eta.^2);             % 距离补偿调频率%%%%%%%%精确的是f_eta为变量
%  计算滤波器
Hsrcf = exp(-1j*pi*f_tau_X.^2.*(1./Ksrc));
%  匹配滤波
Saf_tf_33 =ifft(Saf_tf.*Hsrcf,Nrg,2);

% 绘图：方位向傅里叶变换
figure(2)
subplot(222),imagesc(abs(Saf_tf))
ylabel('方位频率(采样点)'),title('(b)方位向fft后的数据')
subplot(223),imagesc(abs(Saf_tf_33))
xlabel('距离向(采样点)'),ylabel('方位频率(采样点)'),title('(c)SRC后的数据')


%% 信号设置--》距离徙动校正
Drd = sqrt(1-lambda^2*f_etaY.^2/(4*Vr^2));              % 距离多普勒域中的徙动因子
RCM = r_tauX./Drd-r_tauX;                               % 需要校正的距离徙动量
RCM = R0 + RCM - R_eta_c;                               % 将距离徙动量转换到原图像坐标系中
offset = RCM/rho_r;                                     % 将距离徙动量转换为距离单元偏移量
%  计算插值系数表
x_tmp = repmat(-4:3,[16,1]);                            % 插值长度                          
x_tmp = x_tmp + repmat(((1:16)/16).',[1,8]);            % 量化位移
hx = sinc(x_tmp);                                       % 生成插值核
kwin = repmat(kaiser(8,2.5).',[16,1]);                  % 加窗
hx = kwin.*hx;
hx = hx./repmat(sum(hx,2),[1,8]);                       % 核的归一化
%  插值表校正
tic
wait_title = waitbar(0,'开始进行距离徙动校正 ...');  
pause(1);
Srcmf_tf = zeros(Naz,Nrg);
for a_tmp = 1 : Naz
    for r_tmp = 1 : Nrg
        offset_ceil = ceil(offset(a_tmp,r_tmp));
        offset_frac = round((offset_ceil - offset(a_tmp,r_tmp)) * 16);
        if offset_frac == 0
           Srcmf_tf(a_tmp,r_tmp) = Saf_tf_33(a_tmp,ceil(mod(r_tmp+offset_ceil-0.1,Nrg))); 
        else
           Srcmf_tf(a_tmp,r_tmp) = Saf_tf_33(a_tmp,ceil(mod((r_tmp+offset_ceil-4:r_tmp+offset_ceil+3)-0.1,Nrg)))*hx(offset_frac,:).';
        end
    end
    
    pause(0.001);
    Time_Trans   = Time_Transform(toc);
    Time_Disp    = Time_Display(Time_Trans);
    Display_Data = num2str(roundn(a_tmp/Naz*100,-1));
    Display_Str  = ['Computation Progress ... ',Display_Data,'%',' --- ',...
                    'Using Time: ',Time_Disp];
    waitbar(a_tmp/Naz,wait_title,Display_Str)
end
pause(1);
close(wait_title);
toc

%% 信号设置--》方位压缩
%  计算滤波器
Haf = exp(+1j*4*pi*f0*r_tauX.*Drd/c);
Haf_offset = exp(-1j*2*pi*f_etaY.*t_eta_c);
%  匹配滤波
Soutf_tf = Srcmf_tf.*Haf.*Haf_offset;
soutt_tt = ifft(Soutf_tf,Naz,1);
% 绘图：方位压缩
figure(2)
subplot(224),imagesc( abs(soutt_tt))
xlabel('距离向(采样点)'),ylabel('方位向(采样点)'),title('(d)RCMC和方位压缩后的数据')
sgtitle('图6.17 方式2下SRC的精确实现')

%% 信号设置:点目标B的分析
for iiiii=1:1:Ntarget
    len = 32;%点附近len*len框
    cut = -len/2:1:len/2-1;
    n = round(Naz/2 +  NPosition(iiiii,2)/Vr*Fa);
    m = round(Nrg/2 +  2*(NPosition(iiiii,1)-R0)/c*Fr);
    start_tt = soutt_tt(n-len:n+len,m-len:m+len);
    [p,q]=find(abs(start_tt)==max(max(abs(start_tt))));
    start_tt = soutt_tt(n-(len-p)-1 + cut,m-(len-q)-1 + cut);
    
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

%% 时间转换函数
function y = Time_Transform(u)
    Time_in = u(1);
    Hours   = fix(Time_in/3600);
    Minutes = fix((Time_in-Hours*3600)/60);
    Seconds = fix(Time_in-Hours*3600-Minutes*60);
    Time_out = [Hours Minutes Seconds];
    y = Time_out;
end
%% 时间显示函数
function y = Time_Display(u)
    Hours   = u(1);
    Minutes = u(2);
    Seconds = u(3);
    
    if Hours == 0
        if Minutes == 0
            Time_out = [num2str(Seconds),'','s'];
        else
            Time_out = [num2str(Minutes),'','m','',...
                            num2str(Seconds),'','s'];
        end 
    else
        Time_out = [num2str(  Hours),'','h','',...
                        num2str(Minutes),'','m','',...
                        num2str(Seconds),'','s'];
    end
    y = Time_out;
end
