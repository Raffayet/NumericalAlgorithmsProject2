clc;
clear all;

%leva strana diferencijalne jednacine: y'' - 2*y' + y/2*x.^2
%desna strana diferencijalne jednacine: x - 3

levaStrana = @(x, h) NANSLib.ddf(h) - 2*NANSLib.df(h) + x.^2*NANSLib.f(h)/2;
desnaStrana = @(x) x - 3;

%a)

x1 = -1;
x2 = 2;
fX1 = 1;
fX2 = 4;

h = (x2 - x1)/100;

x = x1:h:x2;
fX = NANSLib.finiteDifference(levaStrana, desnaStrana, x1, fX1, x2, fX2, h);
plot(x, fX), hold on

%b)

x = 0.5;
index = floor((x - x1 + (x2 - x1)/100)/(x2 - x1)*100);
scatter(x, fX(index), 'red'), hold on
fX(index)
