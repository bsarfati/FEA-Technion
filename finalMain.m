%Ben Sarfati 885573816
%FEA Final Project
clear; close all; clc

%% General Parameters

%Choice of simulation
problems = {'benchmark','simple','project'}; %Available options
indProblem = 3;

%Choice of element type
elementTypes = {'linear triangular','quadratic triangular'}; %Available options
indElementType = 1;

%Mesh refinement factor (applied to each macro element)
mesh_refinement_factor = 10;

%Conversion from BC label to BC type (1=homog. Dirichlet, 2=homog. Neumann,
%3=prescribed transverse deflection from geometry 1, as below)
labels = {'$w=0$','$\frac{\partial w}{\partial n}=0$','$w=c(-y^2+\frac{1}{4}a^2)$'};

%Plot settings 
fontSize = 30;
axFontSize = 30;
lineWidth = 2;
markerSize = 24;
textSpacing = 0.1;
colors = [viridis(7); viridis(7)];

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
        w = @(x,y) 0; %SOMETHING GOES HERE
end

%% Pre-processor

%Create mesh, assign material properties, and assign loading conditions
switch problem
    case 'benchmark' %(uses modified hw3 geometry)
        %Assing forcing function
        p = @(x,y) p0*ones(size(x));

        %Assign material properties
        k = @(x,y) k0;

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

    case 'simple'   
        %Assing forcing function
        p = @(x,y) f(x,y);

        %Assign material properties
        k = @(x,y) 0;

        %Create node coordinates and macro element
        Vraw = [0 0; 0.5 0; 0.5 0.5; 0.25 0; 0.5 0.25; 0.25 0.25];
        Eraw = 1:6;
        
        %Mesh elements
        [E,N] = meshElement(Vraw(Eraw,:),mesh_refinement_factor);    
    case 'project'
        % [E,N] = [NaN NaN]; %SOMETHING GOES HERE
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
        dirichletBoundaryNodes = boundaryNodes(boundaryLabels==1);
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
            boundaryLabels(indCorners(2)+1:end) = 2;
            boundaryLabels(1:indCorners(1)-1) = 2;
        else
            boundaryLabels(indCorners(2):indCorners(1)) = 2;
            boundaryLabels(indCorners(2)+1:end) = 1;
            boundaryLabels(1:indCorners(1)-1) = 1;
        end
        % % if indCorners(2) > indCorners(1) %switch corners
        % %     boundaryLabels(indCorners(1)+1:indCorners(2)-1) = 1;
        % %     boundaryLabels(indCorners(2):end) = 2;
        % %     boundaryLabels(1:indCorners(1)) = 2;
        % % else
        % %     boundaryLabels(indCorners(2)+1:indCorners(1)-1) = 2;
        % %     boundaryLabels(indCorners(2):end) = 1;
        % %     boundaryLabels(1:indCorners(1)) = 1;
        % % end
        dirichletBoundaryNodes = boundaryNodes(boundaryLabels==1);
    case 'project'
        % SOMETHING GOES HERE
end

%% Visualize mesh

figure; grid; set(gcf,'color','w'); hold on
plotElement(E,N);
labelNums = 1:numBoundaryLabels;
handles = zeros(numBoundaryLabels,1);
for i = 1:numBoundaryNodes
    newLabel = find(boundaryLabels(i) == labelNums);
    if newLabel
        labelNums(newLabel) = 0;
        handles(newLabel) = plot(N(boundaryNodes(i),1),N(boundaryNodes(i),2),'.','markersize',markerSize*2,'Color',colors(boundaryLabels(i)*3,:));
        continue;
    end
    plot(N(boundaryNodes(i),1),N(boundaryNodes(i),2),'.','markersize',markerSize*2,'Color',colors(boundaryLabels(i)*3,:))
end

xlabel('x','FontSize',fontSize,'Interpreter','latex')
ylabel('y','FontSize',fontSize,'Interpreter','latex')
title('Boundary Conditions','FontSize',fontSize)
legend(handles,labels{1:numBoundaryLabels},'Interpreter','latex','FontSize',fontSize,'location','nw')
set(gca,'fontSize',fontSize)

%% Processor

%Retrieve basis functions, gradients, jacobians according to element type
[phi,Bhat] = retrieveMapping(elementTypes{indElementType});

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

    %Calculate local mass matrix using degree 3 Gauss Quadrature
    Me = gaussQuadrature(integrandMe,3);

    %Write integrand for local stiffness matrix
    integrandKe = @(xi_e,eta_e) B(xi_e,eta_e)'*B(xi_e,eta_e)*detJe(xi_e,eta_e);

    %Calculate local stiffness matrix using degree 0 Gauss Quadrature
    Ke = gaussQuadrature(integrandKe,0); %will this break when bhat non-const?

    %Write local forcing function
    fe = p(N(currentE,1),N(currentE,2));

    %"Add" local matrices to global matrices
    M(currentE,currentE) = M(currentE,currentE)+Me;
    K(currentE,currentE) = K(currentE,currentE)+Ke;
    F(currentE) = F(currentE)+Me*fe;
end

%Correct global matrices by adding boundary conditions
K(dirichletBoundaryNodes,:) = 0;
K(:,dirichletBoundaryNodes) = 0; %Optional; results from homog. BC
K(dirichletBoundaryNodes,dirichletBoundaryNodes) = eye(length(dirichletBoundaryNodes));
M(dirichletBoundaryNodes,:) = 0;
M(:,dirichletBoundaryNodes) = 0; %Optional; results from homog. BC
F(dirichletBoundaryNodes) = 0;

%Solve (also possible to omit boundary rows)
switch problem
    case 'benchmark'
        a = (K+k0*M)\F;
    case 'simple'
        a = K\F;
    case 'project'
end

%% Visualize

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

e1 = abs((aTheo-a)./aTheo);
e1(aTheo==0 & a==0) = 0;
lims = max(e1);
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