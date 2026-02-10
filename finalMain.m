%Ben Sarfati 885573816
%FEA Final Project
clear; close all; clc

%% General Parameters

%Choice of simulation
problems = {'benchmark','simple','project'}; %Available options
indProblem = 3;

%Choice of element type
elementTypes = {'linear triangular','quadratic triangular'}; %Available options
indElementType = 2;

%Mesh refinement factor (applied to each macro element)
mesh_refinement_factor = 6;

%Plot settings 
fontSize = 30;
axFontSize = 30;
lineWidth = 2;
markerSize = 15;
textSpacing = 0.1;
colors = orderedcolors("gem");

%Conversion from BC label to BC type (1=homog. Dirichlet, 2=homog. Neumann,
%3,4=prescribed transverse deflection from geometry 1, as below); plot
%labels
labels = {'$w=0$','$\frac{\partial w}{\partial n}=0$','$w=c(-y^2+\frac{1}{4}a^2)$','$w=0.2c(-y^2+\frac{1}{4}a^2)$'};

%% Problem-specific parameters

problem = problems{indProblem};
switch problem
    case 'benchmark'
        %Benchmark problem parameters
        R = 10; %Radius
        k0 = 12.1; %Stiffness
        p0 = 5; %Transverse pressure
    case 'simple'
        f = @(x,y) 2*pi^2*sin(pi*x).*sin(pi*y);
    case 'project'
        Tproj = 1.6;
        kMax = 5;
        a = 3;
        b = 12;
        p0 = 5;
        c = 0.9;
        wBarL = @(y) c*(-y.^2+(a/2)^2);
        wBarR = @(y) 0.2*c*(-y.^2+(a/2)^2);
end

%% Pre-processor

%Select element type 
elementType = elementTypes{indElementType};

%Create mesh, assign material properties, and assign loading conditions
switch problem
    case 'benchmark' %(uses modified hw3 geometry)
        %Assign forcing function
        p = @(x,y) p0*ones(size(x));

        %Assign material properties
        k = @(x,y) k0;

        %Assign leading coefficient
        T = 1;

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
            [Emeshed{i},Nmeshed{i}] = meshElement(Vraw(Eraw{i},:),mesh_refinement_factor,elementType);
        end
        
        %Merge meshes
        E = Emeshed{1};
        N = Nmeshed{1};
        for i = 2:6
            [E,N] = mergeMesh(E,N,Emeshed{i},Nmeshed{i});
        end

    case 'simple'   
        %Assign forcing function
        p = @(x,y) f(x,y);

        %Assign material properties
        k = @(x,y) 0;

        %Assign leading coefficient
        T = 1;

        %Create node coordinates and macro element
        Vraw = [0 0; 0.5 0; 0.5 0.5; 0.25 0; 0.5 0.25; 0.25 0.25];
        Eraw = 1:6;
        
        %Mesh elements
        [E,N] = meshElement(Vraw(Eraw,:),mesh_refinement_factor,elementType);    
    case 'project'
        %Assing forcing function
        p = @(x,y) p0*ones(size(x));

        %Assign material properties
        k = @(x,y) kMax*sin(pi*x/b).*cos(pi*y/a);

        %Assign leading coefficient
        T = Tproj;

        %Create node coordinates and macro element
        Vraw = [0 0; b 0; b a/2; 0 a/2; b/2 0; b a/4; b/2 a/2; 0 a/4; b/2 a/4];
        Eraw = 1:9;

        %Mesh elements
        [E,N] = meshElement(Vraw(Eraw,:),mesh_refinement_factor,elementType);
end

%Retrieve ordered boundary nodes
boundaryNodes = getOrderedBoundary(E);
numBoundaryNodes = length(boundaryNodes);

%Assign boundary labels based on problem BC's
boundaryLabels = zeros(numBoundaryNodes,1);
switch problem
    case 'benchmark' %(uses modified hw3 geometry)
        %Boundary BC is homogeneous
        numBoundaryLabels = 1;
        boundaryLabels = ones(numBoundaryNodes,1);
        type1BoundaryNodes = boundaryNodes(boundaryLabels==1);
        type3BoundaryNodes = [];
        type4BoundaryNodes = [];
    case 'simple'   
        numBoundaryLabels = 2;

        %Locate corners of macro element
        indCorners = zeros(3,1);
        for i = 1:3
            indCorners(i) = find(vecnorm(N(boundaryNodes,:)-Vraw(i,:),2,2) == 0);
        end

        %Use order of occurrence of corners to determine boundary direction
        [~,indNew2IndOld] = sort(indCorners);

        %Make boundary clockwise no matter what
        if indNew2IndOld(mod(find(indNew2IndOld==1),3)+1) == 3
            boundaryNodes = flipud(boundaryNodes);
            indCorners = numBoundaryNodes+1-indCorners;
        end

        %Label; only the bottom BC (incl. corners!) gets Dirichlet
        if indCorners(2) > indCorners(1)
            boundaryLabels(indCorners(1):indCorners(2)) = 1;
            boundaryLabels(indCorners(2)+1:end) = 2; %breaks if indCorners=end
            boundaryLabels(1:indCorners(1)-1) = 2;
        else
            boundaryLabels(indCorners(2)+1:indCorners(1)-1) = 2;
            boundaryLabels(indCorners(1):end) = 1;
            boundaryLabels(1:indCorners(2)) = 1;
        end
        type1BoundaryNodes = boundaryNodes(boundaryLabels==1);
        type3BoundaryNodes = [];
        type4BoundaryNodes = [];
    case 'project'
        numBoundaryLabels = 4;

        %Locate corners of macro element
        indCorners = zeros(4,1);
        for i = 1:4
            indCorners(i) = find(vecnorm(N(boundaryNodes,:)-Vraw(i,:),2,2) == 0);
        end

        %Use order of occurrence of corners to determine boundary direction
        [~,indNew2IndOld] = sort(indCorners);

        %Make boundary clockwise no matter what
        if indNew2IndOld(mod(find(indNew2IndOld==1),3)+1) ~= 2
            boundaryNodes = flipud(boundaryNodes);
            indCorners = numBoundaryNodes+1-indCorners;
        end

        %Label; bottom gets symmetry (lowest prio at corners)
        if indCorners(2) > indCorners(1)
            boundaryLabels(indCorners(1)+1:indCorners(2)-1) = 2;
        else
            boundaryLabels(indCorners(1)+1:end) = 2;
            boundaryLabels(1:indCorners(2)-1) = 2;
        end

        %Right gets wBarR (secondary prio at corners)
        if indCorners(3) > indCorners(2)
            boundaryLabels(indCorners(2):indCorners(3)-1) = 4;
        else
            boundaryLabels(indCorners(2):end) = 4;
            boundaryLabels(1:indCorners(3)-1) = 4;
        end

        %Top gets 0 disp (highest prio. at corners)
        if indCorners(4) > indCorners(3)
            boundaryLabels(indCorners(3):indCorners(4)) = 1;
        else
            boundaryLabels(indCorners(3):end) = 1;
            boundaryLabels(1:indCorners(4)) = 1;
        end

        %Left gets wBarL (secondary prio. at corners)
        boundaryLabels(boundaryLabels==0) = 3;

        type1BoundaryNodes = boundaryNodes(boundaryLabels==1);
        type3BoundaryNodes = boundaryNodes(boundaryLabels==3);
        type4BoundaryNodes = boundaryNodes(boundaryLabels==4);
end

%% Visualize mesh

figure; set(gcf,'color','w'); hold on
plotElement(E,N);
labelNums = 1:numBoundaryLabels;
handles = zeros(numBoundaryLabels,1);
for i = 1:numBoundaryNodes
    newLabel = find(boundaryLabels(i) == labelNums);
    if newLabel
        labelNums(newLabel) = 0;
        handles(newLabel) = plot(N(boundaryNodes(i),1),N(boundaryNodes(i),2),'ko','markersize',markerSize,'markerfaceColor',colors(boundaryLabels(i),:));
        continue;
    end
    plot(N(boundaryNodes(i),1),N(boundaryNodes(i),2),'ko','markersize',markerSize,'markerfaceColor',colors(boundaryLabels(i),:))
end

xlabel('x','FontSize',fontSize,'Interpreter','latex')
ylabel('y','FontSize',fontSize,'Interpreter','latex')
title('Boundary Conditions','FontSize',fontSize)
legend(handles,labels{1:numBoundaryLabels},'Interpreter','latex','FontSize',fontSize,'location','eastoutside')
set(gca,'fontSize',fontSize)

%% Processor

%Write appropriate Gauss quadrature orders according to element type
numGaussPoints = [3 0; 4 3]; %inexact for quad tris+M (4 pts but integrand order 4)

%Retrieve basis functions, gradients, jacobians according to element type
[phi,Bhat] = retrieveMapping(elementType);

%Add contribution of each element to global matrices
M = zeros(size(N,1));
Kk = M;
K = M;
F = zeros(size(N,1),1);
for currentE = E'
    %Calculate Jacobian of current element, and functions of it
    Je = @(xi_e,eta_e) Bhat(xi_e,eta_e)*N(currentE,:);
    detJe = @(xi_e,eta_e) det(Je(xi_e,eta_e));
    B = @(xi_e,eta_e) Je(xi_e,eta_e)\Bhat(xi_e,eta_e);

    %Optional: check jacobians
    for i = 1:length(currentE)
        detJs(i) = detJe(N(currentE(i),1),N(currentE(i),2));
    end
    if any(detJs < 0)
        warning('negative Jacobians detected.');
        plot(N(currentE,1),N(currentE,2),'r*','markersize',markerSize,'HandleVisibility', 'off');
    end

    %Write integrand for local mass matrix
    integrandMe = @(xi_e,eta_e) phi(xi_e,eta_e)'*phi(xi_e,eta_e)*detJe(xi_e,eta_e);

    %Calculate local mass matrix using appropriate degree Gauss Quadrature
    Me = gaussQuadratureTri(integrandMe,numGaussPoints(indElementType,1));

    %Write integrand for local stiffness matrix
    integrandKe = @(xi_e,eta_e) B(xi_e,eta_e)'*B(xi_e,eta_e)*detJe(xi_e,eta_e);

    %Calculate local stiffness matrix using appropriate degree Gauss Quadrature
    Ke = gaussQuadratureTri(integrandKe,numGaussPoints(indElementType,2)); %will this break when bhat non-const?

    %Write local forcing function
    pe = p(N(currentE,1),N(currentE,2));

    %Write mapping to global coords from this element's coords
    x = @(xi_e,eta_e) phi(xi_e,eta_e)*N(currentE,1);
    y = @(xi_e,eta_e) phi(xi_e,eta_e)*N(currentE,2);

    %Write integrand for local massy stiffnessy kappa matrix
    integrandKke = @(xi_e,eta_e) integrandMe(xi_e,eta_e)*k(x(xi_e,eta_e),y(xi_e,eta_e));

    %Calculate local massy stiffnessy kappa matrix using degree 4 Gauss 
    % Quadrature (or higher if available)
    Kke = gaussQuadratureTri(integrandKke,4);

    %"Add" local matrices to global matrices
    Kk(currentE,currentE) = Kk(currentE,currentE)+Kke;
    M(currentE,currentE) = M(currentE,currentE)+Me;
    K(currentE,currentE) = K(currentE,currentE)+Ke;
    F(currentE) = F(currentE)+Me*pe;
end

%Correct global matrices by adding zero type boundary conditions
K(type1BoundaryNodes,:) = 0;
K(:,type1BoundaryNodes) = 0; %Optional; results from homog. BC
K(type1BoundaryNodes,type1BoundaryNodes) = eye(length(type1BoundaryNodes));
M(type1BoundaryNodes,:) = 0;
M(:,type1BoundaryNodes) = 0; %Optional; results from homog. BC
Kk(type1BoundaryNodes,:) = 0;
Kk(:,type1BoundaryNodes) = 0; %Optional; results from homog. BC
F(type1BoundaryNodes) = 0;

%add parabola boundary conditions from project problem
parabolaBoundaryNodes = [type3BoundaryNodes; type4BoundaryNodes];
if ~isempty(parabolaBoundaryNodes)
    K(parabolaBoundaryNodes,:) = 0;
    K(parabolaBoundaryNodes,parabolaBoundaryNodes) = eye(length(parabolaBoundaryNodes));
    M(parabolaBoundaryNodes,:) = 0;
    Kk(parabolaBoundaryNodes,:) = 0;
    F(type3BoundaryNodes) = wBarL(N(type3BoundaryNodes,2));
    F(type4BoundaryNodes) = wBarR(N(type4BoundaryNodes,2)); 
end

%Solve (also possible to omit boundary rows)
a = (T*K+Kk)\F;

%% Visualize (for myself)

switch problem
    case 'benchmark'
        wTheo = @(r) p0/k0*(1-besseli(0,sqrt(k0)*r)/besseli(0,sqrt(k0)*R));
        aTheo = wTheo(vecnorm(N,2,2)); 
    case 'simple'
        wTheo = @(x,y) sin(pi*x).*sin(pi*y);
        aTheo = wTheo(N(:,1),N(:,2)); 
    case 'project'
        %IUoaoweinrfvs
end

if ~strcmp(problem,'project')
    e1 = abs((aTheo-a)./aTheo);
    e1(aTheo==0 & a==0) = 0;
    lims = max(e1);
    figure;
    plot(e1,'*')
    xline(boundaryNodes)
    title('with local forcing vector')
end

%% Post-processing

