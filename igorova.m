clc
clear all

% a)
leva= @(x,h) NANSLib.ddf(h)-2*NANSLib.df(h)+x.^2*NANSLib.f(h)/2;
desna=@(x)x-3;
x1=-1;
fx1=1;
x2=2;
fx2=4;
h=(x2-x1)/1000;
x=x1:h:x2;


fX = NANSLib.finiteDifference(leva,desna,x1,fx1,x2,fx2,h);
plot(x, fX), hold on

% b)
x=0.5;
p=ceil((x-x1)/(x2-x1)*1000);
scatter(x,fX(p),'b'),hold on
fX(p)