clc
clear
close all

% 设置调频斜率K 和 基本时间长度T0
K  = 25;
T0 = 1;

% 图片设置
h         = figure; set( h, 'position', [100,10,600,700]);
sub_row   = 3;
sub_col   = 2;
sub_count = 0;

% 迭代
T = T0;
for times = 1:3
    
    % 计算脉冲宽度
    T = T0 * 2^(times-1);
    % 计算B和TBP
    B   = K*T;
    TBP = B*T;
    
    % 生成抽样时间间隔dt
    fs = 1.25*B;
    dt = 1/fs;
    Nt = ceil(T/dt);
    
    % 对抽样参数进行修正
    Nt = 2^( ceil( log2(Nt) ) );
    dt = T/Nt;
    fs = 1/dt;
    
    % 生成时间序列t和信号st
    t  = -T/2:dt:T/2-dt;
    st = exp( 1i.*pi.*K.*t.^2 );
    
    % 得到频点
    fre = ( 0:Nt-1 )*fs / Nt - fs/2;
    
    % 进行fft 得到频谱sw
    sw = fftshift( fft( fftshift( st ) ) ) * dt;
    
    % 得到驻定相位法得到的频谱sw2
    sw2 = 1./sqrt(abs(K)) .* exp(-1i.*pi.*fre.^2./K) .* ( fre <=B./2 & fre >= -B./2 );
    if K > 0
        sw2 = sw2.*exp(1i.*pi/4);
    else
        sw2 = sw2.*exp(-1i.*pi/4);
    end
    
    % 画图 频谱幅度
    sub_count = sub_count + 1; subplot( sub_row, sub_col, sub_count );
    plot( fre, abs(sw), 'r', fre, abs(sw2), 'b' ); %legend('DFT','驻定相位法');
    ylabel('幅度'); grid on
    if times == 5
       xlabel('频率（Hz）');  
    end
    
    % 画图 频谱相位
    sub_count = sub_count + 1; subplot( sub_row, sub_col, sub_count );
    plot( fre, phase(sw), 'r', fre, phase(sw2), 'b' ); %legend('DFT','驻定相位法');
    ylabel('相位（弧度）'); grid on; title( ['TBP = ', num2str(TBP)] );
    if times == 5
       xlabel('频率（Hz）');  
    end    
    
end
sgtitle('图3.3 不同TBP的DFT频谱变化');

