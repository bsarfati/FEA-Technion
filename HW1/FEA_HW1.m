%% FEA_HW1
%Ben Sarfati 941180069
clear all; close all; clc

%% Q1.2

func = @(x) (x+cos(x)).*exp(x);
n = (0:(1/6):0.5)';
f = func(n);
N = @(x) [(x-n(2)).*(x-n(3)).*(x-n(4))./(n(1)-n(2))./(n(1)-n(3))./(n(1)-n(4))...
    (x-n(1)).*(x-n(3)).*(x-n(4))./(n(2)-n(1))./(n(2)-n(3))./(n(2)-n(4))...
    (x-n(1)).*(x-n(2)).*(x-n(4))./(n(3)-n(1))./(n(3)-n(2))./(n(3)-n(4))...
    (x-n(1)).*(x-n(2)).*(x-n(3))./(n(4)-n(1))./(n(4)-n(2))./(n(4)-n(3))];

xvals = (linspace(0,0.5,10000))';
figure; grid; hold on; set(gcf,'Color',[1 1 1]);
plot(xvals,func(xvals),'b')
plot(xvals,N(xvals)*f,'r')
legend('Analytical','Numerical')
xlabel('x')
ylabel('f(x)')

%% Q2

xvals = linspace(0,1,10000);
figure; grid; hold on; set(gcf,'Color',[1 1 1]);
plot(xvals,erf(xvals),'b')
plot(xvals,MyErf1(xvals),'r')
plot(xvals,MyErf2(xvals),'k')
legend('Analytical','1 Nodal Point','2 Nodal Points','Location','nw')
xlabel('x')
ylabel('erf(x)')

function I = MyErf1(x)
    func = @(t) exp(-t.^2);
    t = 0;
    f = func((t+1)*(x/2));
    w = 2;
    I = (2/sqrt(pi))*(x/2).*(w*f);
end

function I = MyErf2(x)
    func = @(t) exp(-t.^2);
    t = sqrt(1/3)*(-1:2:1)';
    f = func((t+1)*(x/2));
    w = [1 1];
    I = (2/sqrt(pi))*(x/2).*(w*f);
end