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
mesh_refinement_factor = 6;

%Flag for calculating error norms
calcingErrorNorms = true;

%Plot settings 
fontSize = 30;
axFontSize = 30;
lineWidth = 2;
markerSize = 15;
textSpacing = 0.1;
colors = orderedcolors("gem");
defFieldRes = 10; %Deformation field resolution

%Conversion from BC label to BC type (1=homog. Dirichlet, 2=homog. Neumann,
%3,4=prescribed transverse deflection from geometry 1, as below); plot
%labels
labels = {'$w=0$','$\frac{\partial w}{\partial n}=0$','$w=c(-y^2+\frac{1}{4}a^2)$','$w=0.2c(-y^2+\frac{1}{4}a^2)$'};

%% Problem-specific parameters

k0 = 12.1; %Stiffness
p0 = 1; %Transverse pressure
len = 5;
wid = 4;

%% Pre-processor

%Select element type 
elementType = elementTypes{indElementType};

%Assign forcing function
p = @(x,y) p0*ones(size(x));

%Assign material properties
k = @(x,y) k0;

%Assign leading coefficient
T = 1;

%Create node coordinates and macro element
Vraw = [0 0; len 0; len wid; 0 wid; len/2 0; len wid/2; len/2 wid; 0 wid/2; len/2 wid/2];
Eraw = 1:9;

%Mesh elements
[E,N] = meshElement(Vraw(Eraw,:),mesh_refinement_factor,elementType);

%Retrieve ordered boundary nodes
boundaryNodes = getOrderedBoundary(E);
numBoundaryNodes = length(boundaryNodes);

%Boundary BC is homogeneous
numBoundaryLabels = 1;
boundaryLabels = ones(numBoundaryNodes,1);
type1BoundaryNodes = boundaryNodes(boundaryLabels==1);
type3BoundaryNodes = [];
type4BoundaryNodes = [];

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

%Get rid of boundary nodes
Kinner = K;
% % % % % % Kinner(boundaryNodes,:)
% % % 
% % % % %Solve (also possible to omit boundary rows)
% % % % a = (T*K+Kk)\F;
% % % 
% % % %% Visualize (for myself)
% % % 
% % % switch problem
% % %     case 'benchmark'
% % %         aTheo = wTheo(vecnorm(N,2,2)); 
% % %     case 'simple'
% % %         aTheo = wTheo(N(:,1),N(:,2)); 
% % % end
% % % 
% % % if ~strcmp(problem,'project')
% % %     e1 = abs((aTheo-a)./aTheo);
% % %     e1(aTheo==0 & a==0) = 0;
% % %     lims = max(e1);
% % %     figure;
% % %     plot(e1,'*')
% % %     xline(boundaryNodes)
% % %     title('% error across nodes')
% % % end
% % % 
% % % %% Post-processing
% % % 
% % % %Turn off warning that mergeMesh (TA code) triggers to see if negative
% % % %Jacobians exist
% % % warning('off', 'MATLAB:colon:operandsNotRealScalar')
% % % 
% % % %Reflect results across symmetries to generate complete solution
% % % switch problem
% % %     case 'benchmark'
% % %         Ecomplete = E; Ncomplete = N; aComplete = a;
% % %     case 'simple'
% % %         N2 = [N(:,2) N(:,1)];
% % %         [Equarter,Nquarter,aQuarter] = mergeMeshSoln(E,N,a,E,N2,a);
% % %         N3 = [Nquarter(:,1) 1-Nquarter(:,2)];
% % %         [Ehalf,Nhalf,aHalf] = mergeMeshSoln(Equarter,Nquarter,aQuarter,Equarter,N3,aQuarter);
% % %         N4 = [1-Nhalf(:,1) Nhalf(:,2)];
% % %         [Ecomplete,Ncomplete,aComplete] = mergeMeshSoln(Ehalf,Nhalf,aHalf,Ehalf,N4,aHalf);
% % %     case 'project'
% % %         N2 = [N(:,1) -N(:,2)];
% % %         [Ecomplete,Ncomplete,aComplete] = mergeMeshSoln(E,N,a,E,N2,a);
% % % end
% % % numElems = size(Ecomplete,1);
% % % numNodes = size(Ncomplete,1);
% % % 
% % % %Show that interpolation of nodal values of solution is largely sufficient
% % % %for creating a smooth result; but "cheats" around function requirement
% % % figure; set(gcf,'color','w'); trisurf(Ecomplete(:,1:3), Ncomplete(:,1), Ncomplete(:,2), aComplete, 'EdgeColor','none');
% % % view(2)
% % % axis equal
% % % colorbar
% % % title('Contour Plot of Solution Using 1st Nodal Value of Each Element','FontSize',fontSize)
% % % set(gca,'fontSize',fontSize)
% % % figure; set(gcf,'color','w'); trisurf(Ecomplete(:,1:3), Ncomplete(:,1), Ncomplete(:,2), aComplete, 'EdgeColor','none');
% % % view(2)
% % % axis equal
% % % shading interp
% % % colorbar
% % % title('Contour Plot of Solution Using Interpolated Nodal Values','FontSize',fontSize)
% % % set(gca,'fontSize',fontSize)
% % % 
% % % %Mesh each individual element and sample solution at all resulting nodes in
% % % %order to plot a higher resolution deformation field using the actual
% % % %function solution w(x,y)
% % % [Eplot,Nplot,wPlot] = mapLocalSolution(Ncomplete(Ecomplete(1,:),:),defFieldRes,'linear triangular',aComplete(Ecomplete(1,:)));
% % % for i = 2:numElems
% % %     [EplotElem,NplotElem,wPlotElem] = mapLocalSolution(Ncomplete(Ecomplete(i,:),:),defFieldRes,'linear triangular',aComplete(Ecomplete(i,:)));
% % %     [Eplot,Nplot,wPlot] = mergeMeshSoln(Eplot,Nplot,wPlot,EplotElem,NplotElem,wPlotElem);
% % % end
% % % figure; set(gcf,'color','w'); trisurf(Eplot,Nplot(:,1),Nplot(:,2),wPlot,'EdgeColor','none');
% % % view(2)
% % % axis equal
% % % colorbar
% % % title('Contour Plot of Solution Using w(x,y) sampled across elements','FontSize',fontSize)
% % % set(gca,'fontSize',fontSize)
% % % 
% % % %Again but shading for ultimate plot
% % % figure; set(gcf,'color','w'); trisurf(Eplot,Nplot(:,1),Nplot(:,2),wPlot,'EdgeColor','none');
% % % view(2)
% % % axis equal
% % % colorbar
% % % shading interp
% % % title('Contour Plot of Solution Using w(x,y) sampled across elements + Shading Interp','FontSize',fontSize)
% % % set(gca,'fontSize',fontSize)
% % % 
% % % %Find curvature -Lap(w)=(p-kw)/T by evaluating exact functions p,k at
% % % %points where true function w(x,y) was also evaluated
% % % numPlotPts = size(wPlot,1);
% % % pCurv = zeros(numPlotPts,1);
% % % kCurv = zeros(numPlotPts,1);
% % % for i = 1:numPlotPts
% % %     pCurv(i) = p(Nplot(i,1),Nplot(i,2));
% % %     kCurv(i) = k(Nplot(i,1),Nplot(i,2));
% % % end
% % % 
% % % %Plot "ultimate" curvature plot
% % % figure; set(gcf,'color','w'); trisurf(Eplot,Nplot(:,1),Nplot(:,2),(pCurv-kCurv.*wPlot)/T,'EdgeColor','none');
% % % view(2)
% % % axis equal
% % % colorbar
% % % title('Contour Plot of Curvature Using true solution (No Shading Interp)','FontSize',fontSize)
% % % set(gca,'fontSize',fontSize)
% % % 
% % % %% F-load F-reaction question
% % % 
% % % %Retrieve complete boundary nodes
% % % boundaryNodesComplete = getOrderedBoundary(E);
% % % numBoundaryNodesComplete = length(boundaryNodesComplete);
% % % 
% % % %Do quadrature to get F-load and 1st term of F-reaction
% % % Fload = 0;
% % % Freaction1 = 0;
% % % for currentE = E'
% % %     %Calculate Jacobian of current element, and functions of it
% % %     Je = @(xi_e,eta_e) Bhat(xi_e,eta_e)*N(currentE,:);
% % %     detJe = @(xi_e,eta_e) det(Je(xi_e,eta_e));
% % % 
% % %     %Write mapping to global coords from this element's coords
% % %     x = @(xi_e,eta_e) phi(xi_e,eta_e)*N(currentE,1);
% % %     y = @(xi_e,eta_e) phi(xi_e,eta_e)*N(currentE,2);
% % % 
% % %     %Write numerical solution on this element
% % %     w = @(xi_e,eta_e) phi(xi_e,eta_e)*a(currentE);
% % % 
% % %     %Write integrand for p integral
% % %     integrandp = @(xi_e,eta_e) p(x(xi_e,eta_e),y(xi_e,eta_e))*detJe(xi_e,eta_e);
% % % 
% % %     %Calculate integral of p using degree 4 Gauss Quadrature
% % %     Floade = gaussQuadratureTri(integrandp,4);
% % % 
% % %     %Write integrand for kw integral
% % %     integrandkw = @(xi_e,eta_e) k(x(xi_e,eta_e),y(xi_e,eta_e))*w(x(xi_e,eta_e),y(xi_e,eta_e))*detJe(xi_e,eta_e);
% % % 
% % %     %Calculate integral of p using degree 4 Gauss Quadrature
% % %     Freaction1e = gaussQuadratureTri(integrandkw,4);
% % % 
% % %     %Add to totals
% % %     Fload = Fload+Floade;
% % %     Freaction1 = Freaction1+Freaction1e;
% % % end
% % % 
% % % %2nd term of F-reaction 
% % % %???