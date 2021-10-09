clc;
clear all;

%x = sqrt((10 + y'' - 4y*sin(2x))/12)
%x^2 = (10 + y'' - 4y*sin(2x))/12
%12x^2 = 10 + y'' - 4y*sin(2x)
%y'' = 4y*sin(2x) + 12x^2 - 10

diferencijalna = @(x, fX, dfX) 4*fX*sin(2*x) + 12*x.^2 - 10;

%a)

nfX0 = [0, 4];
x1 = -2;
x2 = 0;

h = (x2 - x1)/1000;

x = x1:h:x2;
fX = NANSLib.rk4N(x1, x2, h, nfX0, diferencijalna);
plot(x, fX), hold on

%b)
nfX0 = [0, 5];
x1 = -1;
x2 = 3;             %trojka je stavljena da bi se produzio interval,
                    %da se lepo vidi vrednost y(x) kada je x = 2 
                    

h = (x2 - x1)/1000;

x = x1:h:x2;
fX = NANSLib.rk4N(x1, x2, h, nfX0, diferencijalna);
p = NANSLib.lSquares(x, fX, 7);
x = linspace(x1, x2, 100);
plot(x, polyval(p, x), [x1, x2], [0, 0], "x"), hold on
scatter(2, polyval(p, 2))

polyval(p, 2)
