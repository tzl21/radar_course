function [ img_wk ] = wKA( s0, lambda, Kr, Vr, Fr, PRF, Rref, f_etac,Tr )
%wKA Foucus your SAR data by Omega-K Algorithm
%   
%   s2 the focused image
%   s1 the focused image without stolt map
%   ...
c = 299792458;
f0 = c / lambda;
Fa = PRF;

%% 1. 参考函数相乘
% 构造参考函数
[Naz, Nrg] = size(s0);
f_tau = ifftshift(-Nrg/2:Nrg/2-1) * Fr / Nrg;
f_eta = (ifftshift(-Naz/2:Naz/2-1) * Fa / Naz).';
% 将频率[-Fa/2, Fa/2]映射回其实际（卷绕前）对应的频率
f_eta = f_eta + round((f_etac - f_eta) / Fa) * Fa;
[f_tau_grid, f_eta_grid] = meshgrid(f_tau, f_eta);

clear f_tau; clear f_eta;

theta_ref = 4*pi*Rref / c * sqrt((f0+f_tau_grid).^2 ...
- c^2*f_eta_grid.^2/(4*Vr^2)) + pi*f_tau_grid.^2/Kr;
Href = exp(1j * theta_ref);
clear theta_ref;
S2df = fft2(s0);
clear s0;
% 因为参考函数相乘中有线性相位2*pi*2*Rref/c/D，所以需要将其补偿掉，不然时域就发生了平移
D_fetac_Vref = sqrt(1-c^2*f_etac^2/(4*Vr^2*f0^2));
% 一致参考函数相乘导致的距离时间和方位时间偏移（主要看theta_ref展开时中关于f_tau和f_eta的线性项）
tau_shift = (2*Rref/c/D_fetac_Vref);
eta_shift = Rref * c * f_etac / (2 * Vr^2 * f0 * D_fetac_Vref);
fprintf('tau_shift:%f\n',(tau_shift*Fr));
fprintf('eta_shift:%f\n',eta_shift*Fa);
clear D_fetac_Vref;

Hshift1 = exp(-1j*2*pi*f_tau_grid*tau_shift) ...
    .*exp(1j*2*pi*f_eta_grid*eta_shift);
clear tau_shift; clear eta_shift;
S_RFM = S2df .* Href .* Hshift1;
clear Href; clear Hshift1; clear S2df; 

a_os_r = Fr/abs(Kr*Tr);
N_BW_r = round(Nrg/a_os_r);            % Kr*Tr包含的点数
window_r = ifftshift(kaiser(N_BW_r,2.5)');    % Kaiser窗
window_r = repmat([window_r(1:ceil(N_BW_r/2)),zeros(1,Nrg-N_BW_r),window_r(ceil(N_BW_r/2)+1:N_BW_r)],Naz,1);
S_RFM = S_RFM.*window_r;
clear window_r
% s1 = ifft2(S_RFM);
% figure;
% imagesc(abs(s1)); title('无Stolt插值的压缩目标');

%% 2. Stolt映射
% 计算映射频率偏移量
f_tau1_0 = sqrt((f0 + 0)^2 - c^2*f_eta_grid.^2/(4*Vr^2)) - f0; % 映射后距离向中心频率
% 将频率[-Fr/2, Fr/2]映射回其实际（卷绕前）对应的频率 
f_tau1_grid = f_tau_grid + round((f_tau1_0 - f_tau_grid)/Fr)*Fr;
OFFSET = sqrt((f_tau1_grid + f0).^2 + c^2*f_eta_grid.^2/(4*Vr^2)) - f0 - f_tau_grid;
OFFSET = OFFSET / (Fr/Nrg);
figure();
imagesc(fftshift(OFFSET,2))
clear f_eta_grid; clear f_tau1_grid;

% 下面类似于距离徙动校正进行插值
% 计算插值核系数表
x_tmp = repmat(-4:3, 16, 1);
offset_tmp = (1:16)/16;
x_tmp = x_tmp + repmat(offset_tmp.', 1, 8);
hx = sinc(x_tmp);
x_tmp16 = x_tmp .* 16;
x_tmp16 = round(x_tmp16 + 16 * 8 / 2);
kwin = repmat(kaiser(16*8, 2.5).', 16, 1);
hx = kwin(x_tmp16) .* hx;
hx = hx ./ repmat(sum(hx, 2), 1, 8);

% 驻定相位原理要求Chirp信号s(t)的持续时间为[-T/2, T/2]，尔后我们才能推出书上的
% stolt映射公式。但我们直接按[-T/2, T/2]截取的信号去做DFT，那么得到的信号频谱相比
% s(t)多出了线性相位，因此我们需要将时间轴循环移位，或直接在频域进行线性相位补偿
% 具体参考使用Matlab做FFT变换值得注意的一些问题.md
Hshift2 = exp(-1j*2*pi*(f_tau_grid*(-Nrg/2/Fr)));
clear f_tau_grid;
S_RFM = S_RFM .* Hshift2;

figure();
imagesc(fftshift(abs(S_RFM)));
% 插值映射
Sstolt = zeros(Naz, Nrg);  % 存放Stolt映射后的信号
hwait=waitbar(0,'请等待>>>>>>>>');
for i = 1:Naz
    for j = 1:Nrg
        offset_int = ceil(OFFSET(i,j));
        offset_frac = round((offset_int - OFFSET(i,j)) * 16);
        if offset_frac == 0
            Sstolt(i,j) = S_RFM(i,ceil(mod(j+offset_int-0.1,Nrg)));   % 利用信号数据S1的周期性假定
        else
            Sstolt(i,j) = S_RFM(i, ceil(mod((j+offset_int-4:j+offset_int+3)-0.1,Nrg))) * hx(offset_frac,:).';
        end
        
    end
    if mod(i,200)== 0
    waitbar(i/Naz,hwait,['Stolt映射中: ', num2str(i/Naz*100), '%']);
    end
end
close(hwait);
clear OFFSET;

Sstolt = Sstolt ./ Hshift2;  % 目标点已经被定位在R0,去除前面的移位效果

figure();
imagesc(fftshift(abs(Sstolt)));

clear Hshift2;
img_wk = ifft2(Sstolt);
end

