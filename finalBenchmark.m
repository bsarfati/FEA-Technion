%Ben Sarfati 885573816
%FEA HW2
clear; close all; clc

%% Parameters

%Mesh refinement factor
mesh_refinement_factor = 3;

%Benchmark problem parameters
R = 10; %Radius
k0 = 1; %Stiffness
p0 = 1; %Transverse pressure
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

%% Pre-processor?

%Get boundary nodes; NEEDS WORK KKKKKKKKKK
boundaryEdges = patchBoundary(E);
boundaryNodes = edgeListToCurve(boundaryEdges);

%Retrieve basis functions, gradients, jacobians according to element type
[phi,Bhat] = retrieveMapping('linear triangular');

%Add contribution of each element to global matrices
for i = 1:size(E,1)
    %Calculate Jacobian of current element
    J = Bhat*N(E(i,:),:)????????????;

    %Write integrand for local mass matrix
    integrandMe = @(xi_e,eta_e) phi(xi_e,eta_e)'*phi(xi_e,eta_e)*det(J);
    Me = gaussQuadrature(integrandMe,p);
end

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