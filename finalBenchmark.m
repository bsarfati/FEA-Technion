%Ben Sarfati 885573816
%FEA HW2
clear; close all; clc

%% Parameters

%Mesh refinement factor
mesh_refinement_factor = 10;

%Benchmark problem parameters
R = 10; %Radius
k0 = 12.1; %Stiffness
p0 = 5; %Transverse pressure
% % % % meshBoundarRes = 100; %Number of points on boundary of mesh

%% Mesh benchmark problem (using HW3 geometry)

% % % % thetaVals = linspace(0,2*pi,meshBoundaryRes)

%Create node coordinates for a uniform 3x3 grid centered at 0
[Vx,Vy] = ndgrid(-1:1,-1:1);
Vraw = [Vx(:) Vy(:)];

%Create all necessary node coords from the above node coords
Vraw = [3*Vraw; 5*Vraw([1 3 7 9],:); sqrt(R^2/2)*Vraw([1 2 3 4 6 7 8 9],:); R*Vraw([2 4 6 8],:)];
Eraw = {[1 9 7 5 8 4]... %Left triangle
 [9 1 3 5 2 6]... %Right triangle (correctly wound)
 [7 9 21 19 8 13 25 12 20]... %Top biquad
 [1 7 19 14 4 12 23 10 17]... %Left biquad
 [3 1 14 16 2 10 22 11 15]... %Bottom biquad
 [9 3 16 21 6 11 24 13 18]}; %Right biquad

%Mesh elements
Emeshed = cell(6,1);
Nmeshed = cell(6,1);
for i = 1:6
    [Emeshed{i},Nmeshed{i}] = meshElement(Vraw(Eraw{i},:),mesh_refinement_factor);
end

%Merge meshes
E = Emeshed{1};
N = Nmeshed{1};
for i = 2:6
    [E,N] = mergeMesh(E,N,Emeshed{i},Nmeshed{i});
end

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
f = zeros(size(N,1),1);
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

%Create global forcing vector (alternative)
f = (f+1)*p0; %I think this is wrong, trying to figure out why
% exactly still (why we need to do local forcing function instead)

%Correct global matrices by adding boundary conditions
K(boundaryNodes,:) = 0;
K(:,boundaryNodes) = 0; %Optional; results from homog. BC
K(boundaryNodes,boundaryNodes) = eye(length(boundaryNodes));
M(boundaryNodes,:) = 0;
M(:,boundaryNodes) = 0; %Optional; results from homog. BC
F(boundaryNodes) = 0;
F2 = M*f;
F2(boundaryNodes) = 0; 

%Solve (also possible to omit boundary rows)
a = (K+k0*M)\F;
a2 = (K+k0*M)\(F2); %Not actually sure if constants can go here KKKKKK

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
% ylim([-lims lims])
figure;
plot(e2,'*')
xline(boundaryNodes)
title('with direct global forcing vector of 1s')
% ylim([-lims lims])
figure;
plot((a-a2)./a,'*')
xline(boundaryNodes)
title('diff between')
figure;
plot((M*f-F)./F,'*')
xline(boundaryNodes)
title('diff between the terms themselves')

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