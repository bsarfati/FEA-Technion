%Ben Sarfati 885573816
%FEA HW2
clear; close all; clc

%% Parameters

%Number of elements
ns = [3 10 100]; 

%Plot settings 
fontSize = 30;
axFontSize = 24;
lineWidth = 2;
markerSize = 2;
 
%% Solve and plot

%Constants
E = 200e9;  %Pa
A = 0.15;   %m^2
L = 2;      %m
p = 7850;   %kg/m^3
w = 2*pi;   %rad/sec  

%Plot exact solution
uExact = @(xval) p*w^2*xval.*(L^2-xval.^2/3)/2/E;
% dudxExact = @(xval) rho*w^2*(L^2-xval.^2)/2/E;
xvals = linspace(0,L,1e5)';
figure;
hold on; grid; set(gcf,'color','w')
plot(xvals,uExact(xvals),'k','LineWidth',lineWidth)
xlabel('x','FontSize',fontSize,'Interpreter','latex')
ylabel('u','FontSize',fontSize,'Interpreter','latex')
set(gca,'FontSize',axFontSize)
labels{1} = 'Exact';

%Generate solution for each element count
for n = ns
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
    u = myThomas(K,F,N);

    %Add to plot 
    plot(x,u,'o--','LineWidth',lineWidth,'MarkerSize',markerSize);
    labels{end+1} = ['$n=' num2str(n) '$'];
end

%Add labels
legend(labels,'Interpreter','latex','FontSize',fontSize,'location','nw')

%% Plot L2 and H1 norms

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
