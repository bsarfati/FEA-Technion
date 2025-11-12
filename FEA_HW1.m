%Ben Sarfati 885573816
%FEA HW1
clear; close all; clc

%%
fontSize = 30;
axFontSize = 24;

%% Q1a

f = @(x) sin(x).*(x+cos(x)).^2;
xi = linspace(0,0.5,4)';
fi = f(xi);
Ni = getLagrangePoly(xi);
interpF = @(x) arrayfun(@(xx) Ni(xx)*fi,x);
xLims = [-0.4 1];
xplot = linspace(xLims(1), xLims(2),1e3)';
figure; hold on; grid; set(gcf,'color','w');
plot(xplot,f(xplot),'linewidth',2);
plot(xplot,interpF(xplot),'linewidth',2);
xlim(xLims)
xlabel('$x$','FontSize',fontSize,'Interpreter','latex')
ylabel('$f(x)$','FontSize',fontSize,'Interpreter','latex')
legend('Analytical','Interpolated','FontSize',fontSize,'Interpreter','latex')
ax = gca;
ax.FontSize = axFontSize;

%% Q1b

xLims = [-0.15 0.65];
xplot = linspace(xLims(1), xLims(2),1e3)';
figure; grid on; set(gcf,'color','w');
plot(xplot,abs(f(xplot)-interpF(xplot)),'linewidth',2);
xlim(xLims)
xlabel('$x$','FontSize',fontSize,'Interpreter','latex')
ylabel('$|f(x)-f_{interp}(x)|$','FontSize',fontSize,'Interpreter','latex')
ax = gca;
ax.FontSize = axFontSize;

%% Q2

xvals = linspace(0,1,1e3);
figure; hold on; grid; set(gcf,'Color','w');
plot(xvals,erf(xvals),'linewidth',2)
plot(xvals,MyErf1(xvals),'linewidth',2)
plot(xvals,MyErf2(xvals),'--','linewidth',2)
legend('Analytical','1 Nodal Point','2 Nodal Points','Location','nw','FontSize',fontSize,'Interpreter','latex')
xlabel('$x$','FontSize',fontSize,'Interpreter','latex')
ylabel('erf$(x)$','FontSize',fontSize,'Interpreter','latex')
ax = gca;
ax.FontSize = axFontSize;

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