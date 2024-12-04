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

M=32;
N=32;
theta = pi/24;          % 扭曲或旋转角度
% 旋转信号
z1= zeros(M,N);
for ii = 1:M
    for jj = 1:N
        x = jj-N/2;
        y = (M+1-ii)-M/2;
        xx = round(x*cos(-theta)-y*sin(-theta)+N/2);
        yy = M+1-round(x*sin(-theta)+y*cos(-theta)+M/2);
        if(xx>=1 && xx<= N && yy>=1 && yy<=M)
            z1(ii,jj) = z(yy,xx);
        end
    end
end

len=32;
freq=8;
[row,col]=size(z1);
Data=fft2(z1);
slice_fft=fft2(z1);
fft_zero1=zeros(len,freq*len);
fft_zero2=zeros(freq*len,freq*len);
for i=1:row
     [~,min_loc]=min(abs(slice_fft(i,:)));
     fft_zero1(i,1:min_loc)=slice_fft(i,1:min_loc);
     fft_zero1(i,N*col-(col-min_loc)+1:N*col)=slice_fft(i,min_loc+1:col);
end
for j=1:N*col
     [~,min_loc]=min(abs(fft_zero1(:,j)));
     fft_zero2(1:min_loc,j)=fft_zero1(1:min_loc,j);
     fft_zero2(N*row-(row-min_loc)+1:N*row,j)=fft_zero1(min_loc+1:row,j);
end
dataq=ifft2(fft_zero2);

figure
subplot(1,2,1),imagesc(abs(Data));axis image on
xlabel('水平频率');ylabel('垂直频率');title('（a）幅度谱');
subplot(1,2,2),contour(abs(dataq),15);
xlabel('水平（采样点）');ylabel('垂直（采样点）');title('（b）处理后的点目标');
sgtitle('图2.20 扭曲频谱图中的补零位置');