%% FEA HW2
%Ben Sarfati 941180069

clear all; close all; clc;

%% Q8+9+10

%Choose # of elements
ns = [1 5 10 50 100 500 1000];
% ns = [2];

L = 2; E = 200E9; A = 15E-4; rho = 7850; w = 2*pi;
hs = []; L2s = []; H1s = [];

ufcn = @(xval) rho*w^2*xval.*(L^2-xval.^2/3)/2/E;
dudxfcn = @(xval) rho*w^2*(L^2-xval.^2)/2/E;
xvals = (linspace(0,L,1E5))';
u = ufcn(xvals);
figure(1); set(gcf,'color','w')
plot(xvals,u,'k','linewidth',2)
hold on; grid on
xlabel('x','FontSize',16)
ylabel('u','FontSize',16)
clear labels;
labels{1} = 'Exact';

figure(2); set(gcf,'color','w')
hold on; grid on;

for n = ns

    labels{end+1}= ['n=' num2str(n)];
    
    h = L/n;
    hs(end+1) = h;
    x = (linspace(0,L,n+1))';
    Ke = E*A/h*[1 -1; -1 1];
%     Me = h/2*[2 -1; -1 2];
    Me = h/6*[2 1; 1 2];

    Kg = zeros(n+1);
    Fg = zeros(n+1,1);
    for e = 1:n
        fe = rho*A*w^2*x(e:e+1);
        Fe = Me*fe;
        for i = e:e+1
            for j = e:e+1
                Kg(i,j) = Kg(i,j)+Ke(i+1-e,j+1-e);
            end
            Fg(i) = Fg(i)+Fe(i+1-e);
        end
    end

    Fg(1) = 0;
    Kg(1,:) = eye(1,n+1);
    Kg(:,1) = eye(n+1,1);
    a = MyThomas(Kg,Fg,n+1);
    
    utilde=genutilde(a,xvals,L);
    figure(1);
    plot(xvals,utilde,'linewidth',1.5,'Color',0.9-0.8*length(hs)/length(ns)*[1 1 0])
    
    xvalsL2 = x+h/2;
    xvalsL2(end) = [];
    
    utildeL2 = genutilde(a,xvalsL2,L);
    uL2 = ufcn(xvalsL2);
    L2s(end+1) = sqrt(h*sum((uL2-utildeL2).^2));
    
    xvals1H1 = x+h/2*(1+1/sqrt(3));
    xvals1H1(end) = [];
    xvals2H1 = x+h/2*(1-1/sqrt(3));
    xvals2H1(end) = [];
    utilde1H1 = genutilde(a,xvals1H1,L);
    utilde2H1 = genutilde(a,xvals2H1,L);
    dutildedx1H1 = diff(a)/h;
    dutildedx2H1 = diff(a)/h;
    u1H1 = ufcn(xvals1H1);
    u2H1 = ufcn(xvals2H1);
    dudx1H1 = dudxfcn(xvals1H1);
    dudx2H1 = dudxfcn(xvals2H1);
    H1s(end+1) = sqrt(h/2*sum(...
        (u1H1-utilde1H1).^2+(dudx1H1-dutildedx1H1).^2+...
        (u2H1-utilde2H1).^2+(dudx2H1-dutildedx2H1).^2 ...
        ));
end

figure(1);
legend(labels,'fontsize',16,'location','nw');

if length(ns)>1
    [f1,params1] = fit(log10(hs)',log10(L2s)','poly1');
    [f2,params2] = fit(log10(hs)',log10(H1s)','poly1');
    xvals = linspace(min(log10(hs)),max(log10(hs)));
    fit1 = f1.p1*xvals+f1.p2;
    fit2 = f2.p1*xvals+f2.p2;

    figure(2);  
    plot(log10(hs),log10(L2s),'*r','MarkerSize',12);
    plot(xvals,fit1,'r','LineWidth',2);
    plot(log10(hs),log10(H1s),'*b','MarkerSize',12);
    plot(xvals,fit2,'b','LineWidth',2);
    legend('L_2 Data',['\alpha=' num2str(f1.p1)],'H^1 Data',['\alpha=' num2str(f2.p1)],'fontsize',16,'location','nw')
    xlabel('log_{10}(h)','FontSize',16)
    ylabel('log_{10}(||e||)','FontSize',16)
end

function out = genutilde(a,xvals,L)
    %Outputs the numerical solution represented by a in L evaluated at xvals
    
    n = length(a)-1;
    h = L/n;
    x = (linspace(0,L,n+1))';
    phi = zeros(length(xvals),n+1);
    for idx = 1:length(xvals)
        xval = xvals(idx);
        if xval<x(2) %case i = 1
            phi(idx,1) = 1-(xval-x(1))/h;
        end  
        for i = 2:n
            if x(i-1)<xval & xval<=x(i)
                phi(idx,i) = (xval-x(i-1))/h;
            end
            if x(i)<xval & xval<x(i+1)
                phi(idx,i) = 1-(xval-x(i))/h;
            end    
        end
        if x(n)<xval %case i = n+1
            phi(idx,n+1) = (xval-x(n))/h;
        end
    end
    
    out = phi*a;
end

function out = MyThomas(A,b,n)

%Thomas Algorithm for tridiagonal matrices
A(1,2) = A(1,2)/A(1,1);
b(1) = b(1)/A(1,1);
for i = 2:n-1
    A(i,i+1) = A(i,i+1)/(A(i,i)-A(i,i-1)*A(i-1,i));
    b(i) = (b(i)-A(i,i-1)*b(i-1))/(A(i,i)-A(i,i-1)*A(i-1,i));
end
b(n) = (b(n)-A(n,n-1)*b(n-1))/(A(n,n)-A(n,n-1)*A(n-1,n));

for i=n-1:-1:1
    b(i) = b(i)-A(i,i+1)*b(i+1);
end

out = b;
end
