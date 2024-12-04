clear;close all;clc;
xi = 8:1:15;
N = length(xi);
x = linspace(8, 15, 1000);
K = length(x);
g = [0.2,0.8,0.6,1.5,1.3,-1.0,0.75,1.1];
y = zeros(K);
x0 = 11.7;
for k =1: K
    for i= 1: N
        y(k) = y(k)+g(i)*sinc(x(k)-xi(i));
    end
end
x1 = linspace(-4, 4, 1000);
figure
subplot(3,1,1);
plot(x1,sinc(x1),'b');title('(a)插值核');
subplot(3,1,2);
plot(xi,g,'ro',x,sinc(x-x0),'b',xi,sinc(xi-x0),'*');title('(b)插值运算');
subplot(3,1,3);
plot(xi,g,'ro',x,y);title('(c)插值后的信号');
sgtitle('图2.14 使用sinc函数插值的图');

