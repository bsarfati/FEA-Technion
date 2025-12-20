%% FEA HW3 Q1
%Ben Sarfati 941180069

clear all; close all; clc;

%%

%Choose # of elements
ns = [2 3 6];
% ns = [1 5 10 50 100 500 1000];

%Choose orders of elements (respective to ns array)
rs = [3 2 1];
% rs = 3*ones(size(ns));

num_sims = length(ns);

L = 1; E = 200E9; A0 = 5E-4; k = 1E7; F = 1E6; q = 1E6;
% L = 1; E = 1; A0 = 1; k = 1; F = 1; q = 1;
ufcn = @(xval) 3*q*L^2/(2*k)+F/k-L/(4*A0*E)*(q*xval.^2+2*q*L*xval-(8*L^2*q+4*F)*log(1+(xval/L)));
qfcn = @(xval) -L/(4*A0*E)*(2*q*xval+2*q*L-(8*L^2*q+4*F)*1./((1+(xval/L))*L));
hs = []; L2s = []; H1s = [];

xvals = (linspace(0,L,1E5))';
uvals = ufcn(xvals);
qvals = qfcn(xvals);
figure(1); set(gcf,'color','w')
plot(xvals,uvals,'-k','linewidth',2)
hold on; grid on
xlabel('x','FontSize',16)
ylabel('u','FontSize',16)
figure(2); set(gcf,'color','w')
plot(xvals,qvals,'-k','linewidth',2)
hold on; grid on
xlabel('x','FontSize',16)
ylabel('q','FontSize',16)
clear labels;
labels{1} = 'Exact';

for sim = 1:num_sims
    n = ns(sim);
    r = rs(sim);

    labels{end+1}= ['n=' num2str(n) ', r=' num2str(r)];
    
    N = r*n+1;
    h = L/n;
    hs(end+1) = h;
    x = (linspace(0,L,N))';
    
    Kg = zeros(N);
    Fg = zeros(N,1);
    [phihat,Bhat] = gen_phihat_and_Bhat(r);
    Me_integrand = @(xi_e) h/2*phihat(xi_e)'*phihat(xi_e);
    p = ceil((2*r+1)/2);
    Me = MyGaussQuadrature(Me_integrand,p);
    for e = 1:n
        Ke_integrand = @(xi_e) E*A0/L*Bhat(xi_e)'*Bhat(xi_e).*(xi_e-1+2*n+2*e);
        p = r;
        Ke = MyGaussQuadrature(Ke_integrand,p);
        fe = q*(L+x((r*e+1-r):(r*e+1)));
        Fe = Me*fe;
        Kg((r*e+1-r):(r*e+1),(r*e+1-r):(r*e+1)) = Kg((r*e+1-r):(r*e+1),(r*e+1-r):(r*e+1))+Ke; 
        Fg((r*e+1-r):(r*e+1)) = Fg((r*e+1-r):(r*e+1))+Fe;
    end
       
    %Boundary conditions
    Kg(1) = Kg(1)+k;
    Fg(end) = Fg(end)+F;
    
    switch r
        case 1
            a = MyThomas(Kg,Fg,N);    
        otherwise
            a = Kg\Fg;
    end
    
    [utilde,qtilde] = gen_utilde_and_qtilde(a,xvals,L,r);
    figure(1);
    plot(xvals,utilde,'linewidth',1.5)
    figure(2);
    plot(xvals,qtilde,'linewidth',1.5)

    if rs == rs(1) %we only want to calculate norms for multiple h's/1 r
        x_e = @(xi_e,e) L/(2*n)*(xi_e-1+2*e);
        L2sq = 0;
        H1sq = 0;
        for e = 1:n
            L2sq_integrand = @(xi_e) h/2*(ufcn(x_e(xi_e,e))-phihat(xi_e)*a((r*e+1-r):(r*e+1)))^2;
            L2sq = L2sq+MyGaussQuadrature(L2sq_integrand,4);
            H1sq_integrand = @(xi_e) h/2*(...
                (ufcn(x_e(xi_e,e))-phihat(xi_e)*a((r*e+1-r):(r*e+1)))^2+...
                (qfcn(x_e(xi_e,e))-2/h*Bhat(xi_e)*a((r*e+1-r):(r*e+1)))^2);
            H1sq = H1sq+MyGaussQuadrature(H1sq_integrand,4);
        end
        L2s(end+1) = sqrt(L2sq);
        H1s(end+1) = sqrt(H1sq);
    end
end

figure(1);
legend(labels,'fontsize',16,'location','nw');
figure(2);
legend(labels,'fontsize',16,'location','ne');

if rs == r(1)
    [f1,params1] = fit(log10(hs)',log10(L2s)','poly1');
    [f2,params2] = fit(log10(hs)',log10(H1s)','poly1');
    hvals = linspace(min(log10(hs)),max(log10(hs)));
    fit1 = f1.p1*hvals+f1.p2;
    fit2 = f2.p1*hvals+f2.p2;

    figure(3); set(gcf,'color','w')
    hold on; grid on;
    plot(log10(hs),log10(L2s),'*r','MarkerSize',12);
    plot(hvals,fit1,'r','LineWidth',2);
    plot(log10(hs),log10(H1s),'*b','MarkerSize',12);
    plot(hvals,fit2,'b','LineWidth',2);
    legend('L_2 Data',['\alpha=' num2str(f1.p1)],'H^1 Data',['\alpha=' num2str(f2.p1)],'fontsize',16,'location','nw')
    xlabel('log_{10}(h)','FontSize',16)
    ylabel('log_{10}(||e||)','FontSize',16)
end

%% Functions

function [utilde,qtilde] = gen_utilde_and_qtilde(a,xvals,L,r)
    N = length(a);
    n = (N-1)/r;
    h = L/n;
    x = (linspace(0,L,N))';  
    xi_e = @(e,xval) 2/h*xval-(x(r*e+1-r)+x(r*e+1))/h;
    [phihat,Bhat] = gen_phihat_and_Bhat(r);
    utilde = zeros(size(xvals));
    qtilde = zeros(size(xvals));
            
    idx = 1;
    xval = xvals(idx);
    for e = 1:n-1
        while xval <= x(r*e+1)
            utilde(idx) = phihat(xi_e(e,xval))*a((r*e+1-r):(r*e+1));
            qtilde(idx) = 2/h*Bhat(xi_e(e,xval))*a((r*e+1-r):(r*e+1));
            idx = idx+1;
            xval = xvals(idx);
        end
    end
    %for e = n
    utilde(idx:end) = phihat(xi_e(n,xvals(idx:end)))*a(r*n+1-r:end);
    qtilde(idx:end) = 2/h*Bhat(xi_e(n,xvals(idx:end)))*a(r*n+1-r:end);
end

function [phihat,Bhat] = gen_phihat_and_Bhat(r)
    p = @(xi_e) xi_e.^(0:r);
    dpdx = @(xi_e) [zeros(length(xi_e),1) xi_e.^(0:(r-1))].*(0:r);
    A = p((-1:2/r:1)');
    phihat = @(xi_e) p(xi_e)/A;
    Bhat = @(xi_e) dpdx(xi_e)/A;    
end

function integral = MyGaussQuadrature(integrand,p)
    switch p
        case 1
            u = 0;
            w = 2; 
        case 2
            u = [-sqrt(1/3); sqrt(1/3)];
            w = [1 1];
        case 3
            u = [-sqrt(3/5); 0; sqrt(3/5)];
            w = [5/9 8/9 5/9];    
        case 4
            u = [-sqrt(3/7-2/7*sqrt(5/6)) -sqrt(3/7+2/7*sqrt(5/6)) sqrt(3/7-2/7*sqrt(5/6)) sqrt(3/7+2/7*sqrt(5/6))];
            w = [1/2+sqrt(30)/36 1/2-sqrt(30)/36 1/2+sqrt(30)/36 1/2-sqrt(30)/36];
    end
    
    integral = 0;
    for k = 1:p
        integral = integral+w(k)*integrand(u(k));
    end
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
