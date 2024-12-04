%%  距离多普勒算法（RDA)大斜视角+方式3的二次压缩（精确）
%方式2：在距离频率-方位时域进行SRC，并假定与方位频率基本无关。此SRC滤波器可以合并到距离匹
%配滤波器中，因而是最简单的实现方式。其在参考频率上是精确的，但在其他频率上会存在误差。
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
Delta_f_dop = 80;               % 多普勒带宽
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
B_r = -50; B_a = +50;                                   % B点位置
C_r = +50; C_a = B_a + (C_r-B_r)*tan(theta_r_c);        % C点位置
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
xlabel('距离向(采样点)'),ylabel('方位向(采样点)'),title('(c)相位');
pause(1);
close(wait_title);
toc

%% 信号设置：距离压缩
%  一次距离压缩
window = kaiser(Nrg,2.5)';                              % 时域窗
Window = fftshift(window);                              % 频域窗
%  计算滤波器
Hrf = Window.*exp(+1j*pi*f_tau_X.^2/Kr);  
%  二次距离压缩
Drd0 = sqrt(1-lambda^2*f_eta_c^2/(4*Vr^2));             % 距离多普勒域中的徙动因子
Ksrc = 2*Vr^2*f0^3*Drd0^3/(c*R0*f_eta_c^2);             % 距离补偿调频率
%  计算滤波器
Hsrcf = exp(-1j*pi*f_tau_X.^2/Ksrc);  
%  匹配滤波
Sf_ft = fft(st_tt,Nrg,2);
Srf_ft_33 = Sf_ft.*Hrf.*Hsrcf;
srt_tt_33 = ifft(Srf_ft_33,Nrg,2);
% 绘图：距离压缩后的数据
figure(2)
subplot(221),imagesc(abs(srt_tt_33))
ylabel('方位向(采样点)'),title('(a)距离压缩')






%% 信号设置--》方位向傅里叶变换
Saf_tf = fft(srt_tt_33,Naz,1);
% 绘图：方位向傅里叶变换
figure(2)
subplot(222),imagesc(abs(Saf_tf))
ylabel('方位频率(采样点)'),title('(b)方位FFT后的数据')


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
           Srcmf_tf(a_tmp,r_tmp) = Saf_tf(a_tmp,ceil(mod(r_tmp+offset_ceil-0.1,Nrg))); 
        else
           Srcmf_tf(a_tmp,r_tmp) = Saf_tf(a_tmp,ceil(mod((r_tmp+offset_ceil-4:r_tmp+offset_ceil+3)-0.1,Nrg)))*hx(offset_frac,:).';
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
% 绘图：距离徒动校正
figure(2)
subplot(223),imagesc( abs(Srcmf_tf))
xlabel('距离向(采样点)'),ylabel('方位频率(采样点)'),title('(c)RCMC后的数据')


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
xlabel('距离向(采样点)'),ylabel('方位向(采样点)'),title('(d)方位压缩')
sgtitle('图6.20 方式3下SRC的近似实现')


%% 信号设置:点目标B的分析
len = 64;%点附近len*len框
cut = -len/2:len/2-1;
start_tt = soutt_tt(round(Naz/2 + 1 + NPosition(2,2)/Vr*Fa) + cut,...
                    round(Nrg/2 + 1 + 2*(NPosition(2,1)-R0)/c*Fr) + cut);%NPosition(x,2) x=1为A点目标，2时为B点目标，3时为c点目标
Start_ff = fft2(start_tt);
% 绘图：点目标B的放大图（幅度和频谱）
figure(6) 
subplot(221),imagesc(abs(start_tt))
axis([1 len,1 len])
xlabel('距离向(采样点)'),ylabel('方位向(采样点)'),title('(a)点目标B的幅度')
subplot(222),imagesc(abs(Start_ff))
axis([1 len,1 len])
xlabel('距离频率(采样点)');ylabel('方位频率(采样点)');title('(b)点目标B的频谱');

%对D点目标进行升采样
freq=8;%升采样倍数
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
Z =abs(dataq);
%绘图：升采样后的B点目标
figure(6)
subplot(223),imagesc(abs(dataq))
axis([1 freq*len,1 freq*len])
xlabel('距离时间(采样点)'),ylabel('方位时间(采样点)'),title('(c)升采样后的点目标B')
subplot(224),contour(Z,20);
axis([1 freq*len,1 freq*len])
xlabel('距离时间(采样点)'),ylabel('方位时间(采样点)'),title('(c)升采样后的点目标B等值线')

%距离剖面图和方位剖面图
%下面分别对点目标中心（二维最大值）做行切片，和列切片。每一个切片，都做 16倍 升采样。
%考虑图像扭曲,将目标进行旋转和扭曲，以使大部分旁瓣对齐至水平轴和垂直轴
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%距离切片
%第一步，找出升采样后点目标中心的位置和大小
[row_test,column_test] = size(dataq);
[~,p_test]= max(abs(dataq));
[v_test,~] = max(abs(dataq));
[~,q_test] = max(v_test);
row_test_max = p_test(q_test); %二维矩阵最大值，所在第几行
column_test_max = q_test;      %二维矩阵最大值，所在第几列
s_ac_test_max = v_test(q_test);%矩阵最大值是
%第二步，找出点目标中心左侧2/3长度内的最大值，以此来计算距离旁瓣的扭曲角度
%注意到，如果直接在左侧2/3长度内寻找最大值，有可能会误判到方位旁瓣上的点；
%因此，由于斜视角的关系，左侧的距离旁瓣一定是向上倾斜的；
%因此下面的寻找中，行直接规定为在点目标中心所在行的上方。
c3 = round(2*(column_test_max+1)/3);
[~,p_test_2] = max(abs(dataq(1:row_test_max+1,1:c3)));  
[v_test_2,~] = max(abs(dataq(1:row_test_max+1,1:c3)));
[~,q_test_2] = max(v_test_2);
row_test_max_2 = p_test_2(q_test_2);%二维矩阵点目标中心左侧2/3长度的最大值，所在第几行
column_test_max_2 = q_test_2;%二维矩阵点目标中心左侧2/3长度的最大值，所在第几列
%第三步，计算距离旁瓣扭曲的转角
range_theta = atan( abs((row_test_max_2-row_test_max)/(column_test_max_2-column_test_max)) );
%结果是弧度，下面转成角度
range_theta = range_theta/pi*180;%这是距离旁瓣扭曲的角度
%第四步，将升采样的结果s_ac_test以角range_theta进行旋转（这里是逆时针旋转）
s_ac_range = imrotate(dataq,range_theta,'bilinear');%采用'bilinear'，双线性插值的方式
%第五步，做切片
[~,p_test_range] = max(abs(s_ac_range));    
[v_test_range,~] = max(abs(s_ac_range));  
[~,q_test_range] = max(v_test_range);
row_test_max_range = p_test_range(q_test_range);%旋转后，点目标中心，所在第几行
column_test_max_range = q_test_range; %旋转后，点目标中心，所在第几列
s_ac_test_row_max = s_ac_range(row_test_max_range,1:end);
s_ac_test_row_max_no_abs = abs(s_ac_test_row_max)./max(abs(s_ac_test_row_max));       % 归一化
s_ac_test_row_max_abs = 20*log10(s_ac_test_row_max_no_abs);                % 对数化

%绘图：距离向切片
figure(7)
subplot(221),imagesc(abs(s_ac_range))
axis([1 freq*len,1 freq*len])
xlabel('距离时间(采样点)'),ylabel('方位时间(采样点)'),title('(a)升采样后的点目标B距离向扭转校正')
subplot(222),plot(s_ac_test_row_max_abs)
axis([1 freq*len,-35 0])
xlabel('距离时间(采样点)'),ylabel('方位时间(采样点)'),title('(b)升采样后的点目标B距离向切片幅度')
subplot(223),plot(angle(s_ac_test_row_max))
axis([1 freq*len,-4 4])
xlabel('距离时间(采样点)'),ylabel('方位时间(采样点)'),title('(b)升采样后的点目标B距离向切片相位')

%给出IRW PSLR ISLR
[irw,locleft,locright] = get_irw(s_ac_test_row_max_no_abs);
[pslr] = get_pslr(s_ac_test_row_max_abs);
[islr] = get_islr(s_ac_test_row_max_no_abs,2.25);
fprintf( 'B点距离向IRW PSLR ISLR为 [%+3.3f，%+3.3f dB，%+3.3f dB]\n', irw/freq, pslr,islr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%方位切片
%第一步同上
%第二步，找出点目标中心上方（方位向）2/3长度内的最大值，以此来计算方位旁瓣的扭曲角度
%注意到，如果直接在上方2/3长度内寻找最大值，由于距离旁瓣较强，所以有可能会误判到距离旁瓣上的点；
%因此，由于斜视角的关系，上方的方位旁瓣一定是向右倾斜的；
%因此下面的寻找中，列直接规定为在点目标中心所在列的右侧。
cm3 = round(2*row_test_max/3);
[~,p_test_3] = max(abs(dataq(1:cm3,column_test_max:end))); 
[v_test_3,~] = max(abs(dataq(1:cm3,column_test_max:end)));  
[~,q_test_3] = max(v_test_3);
row_test_max_3 = p_test_3(q_test_3);%二维矩阵点目标中心上方2/3长度的最大值，所在第几行
column_test_max_3 = q_test_3;%二维矩阵点目标中心上方2/3长度的最大值，所在第几列
%第三步，计算方位旁瓣扭曲的转角
azimuth_theta = atan( abs((column_test_max_3)/(row_test_max_3-row_test_max)) );
%结果是弧度，下面转成角度
azimuth_theta = azimuth_theta/pi*180;%这是方位旁瓣扭曲的角度
%第四步，将升采样的结果s_ac_test以角azimitu_theta进行旋转（这里是逆时针旋转）
s_ac_azimuth = imrotate(dataq,azimuth_theta,'bilinear');%采用'bilinear'，双线性插值的方式
%第五步，找出旋转后的最大值中心，并取出相应的列切片
[~,p_test_azimuth] = max(abs(s_ac_azimuth));
[v_test_azimuth,~] = max(abs(s_ac_azimuth));
[~,q_test_azimuth] = max(v_test_azimuth);
row_test_max_azimuth = p_test_azimuth(q_test_azimuth);%旋转后，点目标中心，所在第几行
column_test_max_azimuth = q_test_azimuth;%旋转后，点目标中心，所在第几列。
s_ac_test_column_max = s_ac_azimuth(1:end,column_test_max_azimuth);
s_ac_test_column_max_no_abs=abs(s_ac_test_column_max)./max(abs(s_ac_test_column_max));
s_ac_test_column_max_abs=20*log10(s_ac_test_column_max_no_abs);

%绘图：方位向切片
figure(8)
subplot(221),imagesc(abs(s_ac_azimuth))
axis([1 freq*len,1 freq*len])
xlabel('距离向(采样点)'),ylabel('方位向(采样点)'),title('(a)升采样后的点目标B方位向扭转校正')
subplot(222),plot(s_ac_test_column_max_abs)
axis([1 freq*len,-35 0])
xlabel('距离向(采样点)'),ylabel('方位向(采样点)'),title('(b)升采样后的点目标B方位向切片幅度')
subplot(223),plot(angle(s_ac_test_column_max))
axis([1 freq*len,-4 4])
xlabel('距离向(采样点)'),ylabel('方位向(采样点)'),title('(b)升采样后的点目标B距离向切片相位')

%给出IRW PSLR ISLR
[irw,locleft,locright] = get_irw(s_ac_test_column_max_no_abs);
[pslr] = get_pslr(s_ac_test_column_max_abs);
[islr] = get_islr(s_ac_test_column_max_no_abs,2.25);
fprintf( 'B点方位向IRW PSLR ISLR为 [%+3.3f，%+3.3f dB，%+3.3f dB]\n', irw/freq, pslr,islr);



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