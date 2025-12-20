%% FEA HW4 Q1
%Ben Sarfati 941180069

clear all; close all; clc;

%% P4

%Choose what to solve for
finding_eigenmodes = true; %true for P8

%Choose # of elements
% ns = 4; %for P4-5
% ns = 50; %for P6-7
% ns = 100; %for P8
% ns = [100 50 40 30 20 10 5]; %for P9

% 	Choose numerical scheme parameter
alphas = [0 0.5 1];
alpha = alphas(1);

%Choose time step
deltat = 0.006; %close to t_crit for n=4
% deltat = 1E-5; %close to t_crit for n=50

%Choose number of iterations
num_iterations = 4;

%--------------------------------ALGORITHM---------------------------------
rs = ones(size(ns));
num_sims = length(ns);

L = 1; kappa = @(x) (x<=1/4)*1+((x>1/4).*(x<=1/2))*2+((x>1/2).*(x<=3/4))*3+(x>3/4)*2;
hs = []; evals_by_n = [];

xvals = (linspace(0,L,1E5))';
figure(1); set(gcf,'color','w')
hold on; grid on
xlabel('x','FontSize',16)
ylabel('u','FontSize',16)
% clear labels;
labels = {};

for sim = 1:num_sims
    n = ns(sim);
    r = rs(sim);

    N = r*n+1;
    h = L/n;
    hs(end+1) = h;
    x_e = @(xi_e,e) L/(2*n)*(xi_e-1+2*e);

    Kg = zeros(N);
    Mg = zeros(N);
    Fg = zeros(N,1);
    [phihat,Bhat] = gen_phihat_and_Bhat(r);
    Me_integrand = @(xi_e) h/2*phihat(xi_e)'*phihat(xi_e);
    p = ceil((2*r+1)/2);
    Me = MyGaussQuadrature(Me_integrand,p);
    for e = 1:n
        Ke_integrand = @(xi_e) 2/h*kappa(x_e(xi_e,e))*Bhat(xi_e)'*Bhat(xi_e);
        p = r;
        Ke = MyGaussQuadrature(Ke_integrand,p);
        Kg((r*e+1-r):(r*e+1),(r*e+1-r):(r*e+1)) = Kg((r*e+1-r):(r*e+1),(r*e+1-r):(r*e+1))+Ke; 
        Mg((r*e+1-r):(r*e+1),(r*e+1-r):(r*e+1)) = Mg((r*e+1-r):(r*e+1),(r*e+1-r):(r*e+1))+Me; 
    end

    %Consider boundary conditions
    Minner = Mg(2:end-1,2:end-1);
    Kinner = Kg(2:end-1,2:end-1);
    Finner = Fg(2:end-1)-Kg(2:end-1,[1 end])*[0;1]; %Only applies to time-marching problem

    if finding_eigenmodes
        num_iterations = 3;
        a = zeros(N,N-2);
        [a(2:end-1,:),D] = eig(Kinner\Minner);
        evals = diag(D);
        if n==100
            eval_theo = evals(1);
            for mode = 1:num_iterations
                eval = evals(mode)
                labels{end+1} = ['mode ' num2str(mode)];
                [utilde,~] = gen_utilde_and_qtilde(a(:,mode),xvals,L,r);
                plot(xvals,utilde,'linewidth',2)
            end
        else
            evals_by_n(end+1) = evals(1);
        end
    else
        %Set up a with initial and boundary conditionzs
        a = zeros(N,num_iterations);
        a(N,:) = 1; %includes corner point with intersecting IC/BC; BC overrules

        Kbreve = Minner+alpha*deltat*Kinner;
        Kbar = Minner-deltat*(1-alpha)*Kinner;
        Fbar = deltat*Finner;

        max(abs(eig(Kbreve\Kbar)))
        figure(1);
        for timestep = 2:num_iterations
            labels{end+1} = '';
            % r is always 1
            a(2:end-1,timestep) = MyThomas(Kbreve,Kbar*a(2:end-1,timestep-1)+Fbar,N-2);
            [utilde,~] = gen_utilde_and_qtilde(a(:,timestep),xvals,L,r);
            plot(xvals,utilde,'Color',0.9-0.8*timestep/num_iterations*[1 1 0])
        end
        labels{end} = ['n=' num2str(n) ', \alpha=' num2str(alpha)];
    end
end

if length(ns)>1
    hs(1) = [];
    eval_errors = (evals_by_n-eval_theo)/eval_theo;
    [f1,params1] = fit(log10(hs)',log10(eval_errors)','poly1');
    hvals = linspace(min(log10(hs)),max(log10(hs)));
    fit1 = f1.p1*hvals+f1.p2;

    figure(3); set(gcf,'color','w')
    hold on; grid on;
    plot(log10(hs),log10(eval_errors),'*r','MarkerSize',12);
    plot(hvals,fit1,'r','LineWidth',2);
    legend('1st Eigenvalue Data',['\alpha=' num2str(f1.p1)],'fontsize',16,'location','nw')
    xlabel('log_{10}(h)','FontSize',16)
    ylabel('log_{10}(||e||)','FontSize',16)
end

%% P5

if ~finding_eigenmodes
    n = 4;
    r = 1;
    N = r*n+1;
    h = L/n;
    x_e = @(xi_e,e) L/(2*n)*(xi_e-1+2*e);

    Kg = zeros(N);
    Fg = zeros(N,1);
    [phihat,Bhat] = gen_phihat_and_Bhat(r);
    for e = 1:n
        Ke_integrand = @(xi_e) 2/h*kappa(x_e(xi_e,e))*Bhat(xi_e)'*Bhat(xi_e);
        p = r;
        Ke = MyGaussQuadrature(Ke_integrand,p);
        Kg((r*e+1-r):(r*e+1),(r*e+1-r):(r*e+1)) = Kg((r*e+1-r):(r*e+1),(r*e+1-r):(r*e+1))+Ke; 
    end

    %Consider boundary conditions
    Kinner = Kg(2:end-1,2:end-1);
    Finner = Fg(2:end-1)-Kg(2:end-1,[1 end])*[0;1];

    aSS = [0; Kinner\Finner; 1];
    [utilde,qtilde] = gen_utilde_and_qtilde(aSS,xvals,L,r);
    plot(xvals,utilde,'r')
    labels{end+1}= 'steady state';
end

figure(1);
legend(labels,'fontsize',16,'location','nw');

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
