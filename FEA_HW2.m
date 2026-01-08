%Ben Sarfati 885573816
%FEA HW2
clear; close all; clc

%% Parameters

%Number of elements
% ns = [3 10 100]; 
ns = [1 5 10 50 100 500 1000 2000]; %For part 3

%Plot settings 
fontSize = 30;
axFontSize = 24;
lineWidth = 2;
markerSize = 24;
 
%% Solve and plot

%Constants
E = 200e9;  %Pa
A = 0.15;   %m^2
L = 2;      %m
p = 7850;   %kg/m^3
w = 2*pi;   %rad/sec  

%Plot exact solution
num = length(ns);
uExact = @(xVal) p*w^2*xVal.*(L^2-xVal.^2/3)/2/E;
xExact = linspace(0,L,1e5)';
figure(1);
hold on; grid; set(gcf,'color','w')
plot(xExact,uExact(xExact),'k','LineWidth',lineWidth)
xlabel('x','FontSize',fontSize,'Interpreter','latex')
ylabel('u','FontSize',fontSize,'Interpreter','latex')
set(gca,'FontSize',axFontSize)
labels{1} = 'Exact';

%Bonus: plot L2 and H1 norms
dudxExact = @(xVal) p*w^2*(L^2-xVal.^2)/2/E;
hs = zeros(num,1); 
L2s = zeros(num,1); 
H1s = zeros(num,1);
figure(2); 
hold on; grid; set(gcf,'color','w')

%Generate solution for each element count
for i = 1:num
    n = ns(i);            %Number of elements
    N = n+1;              %Number of nodes
    x = linspace(0,L,N)'; %Positions of uniformly spaced nodes
    h = L/n;              %Element length
    J = h/2;              %Jacobian for change of var. from glob. to loc.
    B1 = -1/2;            %Local derivative of local shape function 1
    B2 = 1/2;             %Local derivative of local shape function 2
    K = zeros(N);         %Global stiffness matrix initialized
    F = zeros(N,1);       %Global forcing vector initialized
    Ke = E*A/h*[1 -1      %Local stiffness matrix
                -1 1];    
    Me = h/6*[2 1         %Local mass matrix
              1 2];     

    %Generate global stiffness matrix and forcing vector
    for e = 1:n
        %Define local forcing function
        fe = p*A*w^2*[x(e);x(e+1)];
    
        %Get local forcing vector
        Fe = Me*fe;
        
        %Fill global matrices
        K(e:e+1,e:e+1) = K(e:e+1,e:e+1)+Ke;
        F(e:e+1) = F(e:e+1)+Fe;
    end
    
    %Correct K to satisfy homogeneous dirichlet boundary condition
    K(1,:) = eye(1,N);
    F(1) = 0;
    
    %Simplify K using homogeneous dirichlet boundary condition (optional)
    K(:,1) = eye(N,1);
    
    %Solve (also possible using K\F, and also possible to skip a_1)
    a = myThomas(K,F,N);

    %Add to plot 
    figure(1)
    plot(x,a,'--','LineWidth',lineWidth);
    labels{end+1} = ['$n=' num2str(n) '$'];

    %Generate lists of order 1 and 2 Gauss points for each element
    xL2 = x(1:end-1)+h/2;
    x1H1 = x(1:end-1)+h/2*(1+1/sqrt(3));
    x2H1 = x(1:end-1)+h/2*(1-1/sqrt(3));

    %Evaluate the exact solution at Gauss points of every element...
    uExactL2 = uExact(xL2);
    uExact1H1 = uExact(x1H1);
    uExact2H1 = uExact(x2H1);
    
    %...and exact derivative (at 2nd order Gauss points only)
    dudxExact1H1 = dudxExact(x1H1);
    dudxExact2H1 = dudxExact(x2H1);

    %Evaluate the approximate sol'n at midpoint of every element...
    uTildeL2 = (a(1:end-1)+a(2:end))/2;

    %...at 2nd order Gauss points using linear interpolation...
    uTilde1H1 = (1-1/sqrt(3))/2*a(1:end-1)+(1+1/sqrt(3))/2*a(2:end);
    uTilde2H1 = (1+1/sqrt(3))/2*a(1:end-1)+(1-1/sqrt(3))/2*a(2:end);

    %...and evaluate derivative at 2nd order points (piecewise constant)
    dudxTilde1H1 = diff(a)/h;
    dudxTilde2H1 = diff(a)/h;
    
    %Bonus: L2 and H2 norm calculations using Gauss quadrature
    hs(i) = h;
    L2s(i) = sqrt(h*sum((uExactL2-uTildeL2).^2));
    H1s(i) = sqrt(h/2*sum(...
    (uExact1H1-uTilde1H1).^2+(dudxExact1H1-dudxTilde1H1).^2+...
    (uExact2H1-uTilde2H1).^2+(dudxExact2H1-dudxTilde2H1).^2));
end

%Add labels
figure(1)
legend(labels,'Interpreter','latex','FontSize',fontSize,'location','nw')

%Bonus: L2 and H2 norm plots
if length(ns)>1
    [f1,params1] = fit(log10(hs),log10(L2s),'poly1');
    [f2,params2] = fit(log10(hs),log10(H1s),'poly1');
    xValsNorms = linspace(min(log10(hs)),max(log10(hs)));
    fit1 = f1.p1*xValsNorms+f1.p2;
    fit2 = f2.p1*xValsNorms+f2.p2;

    figure(2);  
    plot(log10(hs),log10(L2s),'*r','MarkerSize',markerSize);
    plot(xValsNorms,fit1,'r','LineWidth',2);
    plot(log10(hs),log10(H1s),'*b','MarkerSize',markerSize);
    plot(xValsNorms,fit2,'b','LineWidth',2);
    legend('$L_2$ Data',['$\alpha=' num2str(f1.p1) '$'],'$H^1$ Data',['$\alpha=' num2str(f2.p1) '$'],...
        'Interpreter','latex','fontsize',fontSize,'location','nw')
    xlabel('$\log_{10}(h)$','FontSize',fontSize,'Interpreter','latex')
    ylabel('$\log_{10}(||e||)$','FontSize',fontSize,'Interpreter','latex')
end

%% Functions

function out = myThomas(A,b,n)
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
