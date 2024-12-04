clear;clc;close all;
x=linspace(-13,13,32);
y=linspace(-13,13,32);
for i=1:length(x)
    for j=1:length(y)      
        x1=x(i)*pi;
        y1=y(j)*pi;
        x2=sin(x1)/x1;
         y2=sin(y1)/y1;
        z(i,j)=x2*y2;
    end
end

len=32;
freq=8;
Data=fft2(z);
S_ac_test_2=fft2(z);
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



x = 1:1:freq*len;
y = 1:1:freq*len;
[X,Y] = meshgrid(x,y);   
Z =abs(dataq);
Z(find(abs(dataq)<0.001))=0;
figure
subplot(1,2,1),imagesc(abs(Data));axis image on
xlabel('水平频率');ylabel('垂直频率');title('（a）幅度谱');
subplot(1,2,2),contour(X,Y,Z,20);
xlabel('水平（采样点）');ylabel('垂直（采样点）');title('（b）处理后的点目标');
sgtitle('图2.18 二维点目标分析图示');

dataq_a=dataq(:,125);
dataq_a_max = abs(dataq_a)./max(abs(dataq_a));
dataq_a_max_abs=20*log10(dataq_a_max);

dataq_r=dataq(125,:);
dataq_r_max = abs(dataq_r)./max(abs(dataq_r));
dataq_r_max_abs=20*log10(dataq_r_max);

figure
subplot(2,2,1),plot(dataq_a_max_abs);ylim([-40 0]);
xlabel('水平采样点');ylabel('幅度（dB）');title('（a）水平切片');
subplot(2,2,2),plot(dataq_r_max_abs);ylim([-40 0]);
xlabel('垂直采样点');ylabel('幅度（dB）');title('（b）垂直切片');
subplot(2,2,3),plot(angle(dataq_a));ylim([-4 4]);
xlabel('水平采样点');ylabel('相位（rad）');title('（c）水平切片');
subplot(2,2,4),plot(angle(dataq_r));ylim([-4 4]);
xlabel('垂直采样点');ylabel('相位（rad）');title('（d）垂直切片');
sgtitle('图2.19 点目标一维剖面图');


