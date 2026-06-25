%Ben Sarfati 885573816
%FEA Final Project
clear; close all; clc

%% General Parameters

%Choice of simulation
problem = 'eigenvalue'; %Available options

%Choice of element type
elementTypes = {'linear triangular','quadratic triangular'}; %Available options
indElementType = 1;

%Mesh refinement factor (applied to each macro element)
mesh_refinement_factors = 4:2:24;

%Number of eigenmodes to plot (max 21)
numModes = 4;

%Plot settings 
fontSize = 30;
axFontSize = 30;
lineWidth = 2;
markerSize = 15;
textSpacing = 0.1;
colors = [orderedcolors("gem"); orderedcolors("gem"); orderedcolors("gem")];
defFieldRes = 10; %Deformation field resolution

%Conversion from BC label to BC type (1=homog. Dirichlet, 2=homog. Neumann,
%3,4=prescribed transverse deflection from geometry 1, as below); plot
%labels
labels = {'$w=0$','$\frac{\partial w}{\partial n}=0$','$w=c(-y^2+\frac{1}{4}a^2)$','$w=0.2c(-y^2+\frac{1}{4}a^2)$'};

%% Problem-specific parameters

c = 2.6;
Lx = 5;
Ly = 4;

%% Pre-processor

%Select element type 
elementType = elementTypes{indElementType};

%Assign leading coefficient
T = -c^2;

allEvals = [];
hs = [];
%%%
for mesh_refinement_factor = mesh_refinement_factors
%%%
disp(['mrf' num2str(mesh_refinement_factor)])

%Create node coordinates and macro element
Vraw = [0 0; Lx 0; Lx Ly; 0 Ly; Lx/2 0; Lx Ly/2; Lx/2 Ly; 0 Ly/2; Lx/2 Ly/2];
Eraw = 1:9;

%Mesh elements
[E,N] = meshElement(Vraw(Eraw,:),mesh_refinement_factor,elementType);

%Retrieve ordered boundary nodes
boundaryNodes = getOrderedBoundary(E);
numBoundaryNodes = length(boundaryNodes);

%Retrieve inner nodes 
innerNodes = setdiff(1:size(N,1),boundaryNodes);

%Boundary BC is homogeneous
numBoundaryLabels = 1;
boundaryLabels = ones(numBoundaryNodes,1);

%% Visualize mesh

% figure; set(gcf,'color','w'); hold on
% plotElement(E,N);
% labelNums = 1:numBoundaryLabels;
% handles = zeros(numBoundaryLabels,1);
% for i = 1:numBoundaryNodes
%     newLabel = find(boundaryLabels(i) == labelNums);
%     if newLabel
%         labelNums(newLabel) = 0;
%         handles(newLabel) = plot(N(boundaryNodes(i),1),N(boundaryNodes(i),2),'ko','markersize',markerSize,'markerfaceColor',colors(boundaryLabels(i),:));
%         continue;
%     end
%     plot(N(boundaryNodes(i),1),N(boundaryNodes(i),2),'ko','markersize',markerSize,'markerfaceColor',colors(boundaryLabels(i),:))
% end
% 
% xlabel('x','FontSize',fontSize,'Interpreter','latex')
% ylabel('y','FontSize',fontSize,'Interpreter','latex')
% title('Boundary Conditions','FontSize',fontSize)
% legend(handles,labels{1:numBoundaryLabels},'Interpreter','latex','FontSize',fontSize,'location','eastoutside')
% set(gca,'fontSize',fontSize)

%% Processor

%Retrieve basis functions, gradients, jacobians according to element type
[phi,Bhat] = retrieveMapping(elementType);

%Add contribution of each element to global matrices
M = zeros(size(N,1));
K = M;
areaSum = 0;
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

    %Calculate local mass matrix using high degree Gauss Quadrature
    Me = gaussQuadratureTri(integrandMe,4);

    %Write integrand for local stiffness matrix
    integrandKe = @(xi_e,eta_e) B(xi_e,eta_e)'*B(xi_e,eta_e)*detJe(xi_e,eta_e);

    %Calculate local stiffness matrix using high degree Gauss Quadrature
    Ke = gaussQuadratureTri(integrandKe,4); %will this break when bhat non-const?

    %"Add" local matrices to global matrices
    M(currentE,currentE) = M(currentE,currentE)+Me;
    K(currentE,currentE) = K(currentE,currentE)+Ke;

    %Calculate area via quadrature of determinant
    areaSum = areaSum+gaussQuadratureTri(detJe,4);
end

%Get rid of boundary nodes
Kinner = K;
Minner = M;
Kinner(boundaryNodes,:) = [];
Kinner(:,boundaryNodes) = [];
Minner(boundaryNodes,:) = [];
Minner(:,boundaryNodes) = [];

%Solve eigenvalue problem
[V,D] = eig(Minner, T*Kinner);
[evals,indEvals] = sort(diag(D),'ascend');
evecs = V(:,indEvals(1:numModes));

%Rebuild nodal solution for first four eigenmodes
emodes = zeros(size(N,1),numModes);
emodes(innerNodes,:) = evecs(:,1:numModes);

%Add these evals to the big list
allEvals = [allEvals -evals(1:numModes)]; %idk why the evals are neg
%Calculate characteristic element size; sqrt(average area of elements)
hs(end+1) = sqrt(areaSum/size(E,1));

%%%
end
%%%

%% Post-processing

%Mesh each individual element and sample solution at all resulting nodes in
%order to plot a higher resolution deformation field using the actual
%function solution w(x,y)
[Eplot,Nplot,wPlot] = mapLocalSolution(N(E(1,:),:),defFieldRes,'linear triangular',emodes(E(1,:),:));
for i = 2:size(E,1)
    [EplotElem,NplotElem,wPlotElem] = mapLocalSolution(N(E(i,:),:),defFieldRes,'linear triangular',emodes(E(i,:),:));
    [Eplot,Nplot,wPlot] = mergeMeshSoln(Eplot,Nplot,wPlot,EplotElem,NplotElem,wPlotElem);
end

for ev = 1:numModes
    figure; set(gcf,'color','w'); trisurf(Eplot,Nplot(:,1),Nplot(:,2),wPlot(:,ev),'EdgeColor','none');
    view(2)
    axis equal
    colorbar
    % shading interp
    title(['Eigenmode ' num2str(ev)],'FontSize',fontSize)
    set(gca,'fontSize',fontSize)
end

%Convergence log-log graph
figure; 
labels = {};
hold on; grid; set(gcf,'color','w')
for i = 1:numModes
    plot(log10(hs),log10(allEvals(i,:)),'o--','MarkerSize',markerSize,'MarkerFaceColor',colors(i,:));
    labels{end+1} = ['$\lambda_' num2str(i) '$'];
end
legend(labels,...
    'Interpreter','latex','fontsize',fontSize,'location','nw')
xlabel('$\log_{10}(h)$','FontSize',fontSize,'Interpreter','latex')
ylabel('$\log_{10}(\lambda)$','FontSize',fontSize,'Interpreter','latex')
set(gca,'fontSize',fontSize)

%Same thing but for frequencies
allFreqs = 1./sqrt(allEvals);
figure; 
labels = {};
hold on; grid; set(gcf,'color','w')
for i = 1:numModes
    plot(log10(hs),log10(allFreqs(i,:)),'o--','MarkerSize',markerSize,'MarkerFaceColor', colors(i,:));
    labels{end+1} = ['$\omega_' num2str(i) '$'];
end
legend(labels,...
    'Interpreter','latex','fontsize',fontSize,'location','nw')
xlabel('$\log_{10}(h)$','FontSize',fontSize,'Interpreter','latex')
ylabel('$\log_{10}(\omega)$','FontSize',fontSize,'Interpreter','latex')
set(gca,'fontSize',fontSize)

% Look for theoretical values within our values
fRatios = zeros(size(allFreqs));
for i = 1:size(allFreqs,2)
    fRatios(:,i) = allFreqs(:,i)/allFreqs(1,i);
end
fRatios