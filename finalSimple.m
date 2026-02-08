%Ben Sarfati 885573816
%FEA HW2
clear; close all; clc

%% Parameters

%Mesh refinement factor
mesh_refinement_factor = 10;

%% Mesh simple problem

%Create node coordinates
Vraw = [0 0; 0.5 0; 0.5 0.5; 0.25 0; 0.5 0.25; 0.25 0.25];

%Create all necessary node coords from the above node coords
Eraw = 1:6;

%Mesh elements
[E,N] = meshElement(Vraw(Eraw,:),mesh_refinement_factor);

%% Vis

figure; hold on;
plotElement(E,N);

%% Pre-processor? Processor?

%Get boundary nodes; NEEDS WORK KKKKKKKKKK
boundaryEdges = patchBoundary(E);
boundaryNodes = edgeListToCurve(boundaryEdges);

%Retrieve basis functions, gradients, jacobians according to element type
[phi,Bhat] = retrieveMapping('linear triangular');

%Add contribution of each element to global matrices
M = zeros(size(N,1));
K = M;
F = zeros(size(N,1),1);
for currentE = E'
    %Calculate Jacobian of current element, and functions of it
    Je = @(xi_e,eta_e) Bhat(xi_e,eta_e)*N(currentE,:);
    detJe = @(xi_e,eta_e) det(Je(xi_e,eta_e));
    B = @(xi_e,eta_e) Je(xi_e,eta_e)\Bhat(xi_e,eta_e);

    %Write integrand for local mass matrix
    integrandMe = @(xi_e,eta_e) phi(xi_e,eta_e)'*phi(xi_e,eta_e)*detJe(xi_e,eta_e);

    %Calculate local mass matrix using degree 2 Gauss Quadrature
    Me = gaussQuadrature(integrandMe,2);

    %Write integrand for local stiffness matrix
    integrandKe = @(xi_e,eta_e) B(xi_e,eta_e)'*B(xi_e,eta_e)*detJe(xi_e,eta_e);

    %Calculate local stiffness matrix using degree 0 Gauss Quadrature
    % Ke = gaussQuadrature(integrandKe,0); %will this break when bhat non-const?
    Ke = gaussQuadrature(integrandKe,1); %will this break when bhat non-const?

    %Write local forcing function
    fe = p0*ones(length(currentE),1);

    %"Add" local matrices to global matrices
    M(currentE,currentE) = M(currentE,currentE)+Me;
    K(currentE,currentE) = K(currentE,currentE)+Ke;
    F(currentE) = F(currentE)+Me*fe;
end

%Correct global matrices by adding boundary conditions
K(boundaryNodes,:) = 0;
K(:,boundaryNodes) = 0; %Optional; results from homog. BC
K(boundaryNodes,boundaryNodes) = eye(length(boundaryNodes));
M(boundaryNodes,:) = 0;
M(:,boundaryNodes) = 0; %Optional; results from homog. BC
F(boundaryNodes) = 0;

%Solve (also possible to omit boundary rows)
a = (K+k0*M)\F;

%% Visualize

wTheo = @(r) p0/k0*(1-besseli(0,sqrt(k0)*r)/besseli(0,sqrt(k0)*R));
aTheo = wTheo(vecnorm(N,2,2)); 

e1 = abs(aTheo-a)./aTheo;
e2 = abs(aTheo-a2)./aTheo;
close all;
lims = max((aTheo-a)./aTheo);
figure;
plot(e1,'*')
xline(boundaryNodes)
title('with local forcing vector')

%% Functions 
function plotElement(elements, nodes)
        scatter(nodes(:,1),nodes(:,2), 'ok');
        t = 0.8;
        for ii=1:size(elements,1)
            centroid = mean(nodes(elements(ii,:),:),1);
            text(centroid(1),centroid(2),num2str(ii),'Color','r');
            for jj=1:size(elements,2)
                node_pos = nodes(elements(ii,jj),:);
                text_pos = t*node_pos+(1-t)*centroid;
                text(text_pos(1),text_pos(2),num2str(elements(ii,jj)),'Color','b');
            end
            pos_array = nodes(elements(ii,:),:);
            pos_array(end+1,:) = pos_array(1,:); %close loop
            plot(pos_array(:,1),pos_array(:,2),'-k')
            
        end
        
end

function integral = gaussQuadrature(integrand,p)
    switch p
        case 0 %Not actual quadrature; just constant val*area
            integral = integrand(NaN,NaN)/2;
            return
        case 1 %might be sufficiently low error order introduced 
            u = [1/3 1/3];
            w = 1;
        case 2
            u = [1/6 1/6; 2/3 1/6; 1/6 2/3];
            w = [1/3 1/3 1/3];
    end
    % switch p
    %     case 1
    %         u = 0;
    %         w = 2; 
    %     case 2
    %         u = [-sqrt(1/3); sqrt(1/3)];
    %         w = [1 1];
    %     case 3
    %         u = [-sqrt(3/5); 0; sqrt(3/5)];
    %         w = [5/9 8/9 5/9];    
    %     case 4
    %         u = [-sqrt(3/7-2/7*sqrt(5/6)) -sqrt(3/7+2/7*sqrt(5/6)) sqrt(3/7-2/7*sqrt(5/6)) sqrt(3/7+2/7*sqrt(5/6))];
    %         w = [1/2+sqrt(30)/36 1/2-sqrt(30)/36 1/2+sqrt(30)/36 1/2-sqrt(30)/36];
    % end
    
    integral = 0;
    for k = 1:p
        integral = integral+w(k)*integrand(u(k,1),u(k,2));
    end
end