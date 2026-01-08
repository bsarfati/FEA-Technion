%Ben Sarfati 885573816
%FEA HW3
clear; close all; clc

%% Parameters 

%Plot settings 
fontSize = 30;
axFontSize = 30;
lineWidth = 2;
markerSize = 24;
textSpacing = 0.1;

%% 1.1

%Create node coordinates for a uniform 3x3 grid centered at 0
[Vx,Vy] = ndgrid(-1:1,-1:1);
V = [Vx(:) Vy(:)];

%Create biquadratic macro element with right-hand node-numbering
E = [1 3 9 7 2 6 8 4 5];

%Plot in the style of plotElement but with local node numbering
plot(V(:,1),V(:,2), '.k','MarkerSize',markerSize);
set(gcf,'color','w')
for i = 1:9
    text(V(E(i),1)+textSpacing,V(E(i),2)+textSpacing,num2str(i),'Color','b','FontSize',fontSize);
end
xlabel('$\xi$','FontSize',fontSize,'Interpreter','latex')
ylabel('$\eta$','FontSize',fontSize,'Interpreter','latex')
set(gca,'FontSize',axFontSize)
xlim([-2 2])
ylim([-2 2])


%% 2.1

%Create all necessary node coords from the above node coords
V = [3*V; 5*V([1 3 7 9],:); 7*V([1 2 3 4 6 7 8 9],:); 10*V([2 4 6 8],:)];
plotNodes(V);
E = {[1 9 7 5 8 4]... %Left triangle
     [1 9 3 5 6 2]... %Right triangle (wound clockwise)
     [7 9 21 19 8 13 25 12 20]... %Top biquad
     [1 7 19 14 4 12 23 10 17]... %Left biquad
     [3 1 14 16 2 10 22 11 15]... %Bottom biquad
     [9 3 16 21 6 11 24 13 18]}; %Right biquad
% plotMesh(E,V);

% testMacro = V(E{3},:);
macro1 = V(E{1},:);
macro3 = V(E{3},:);
[Enew1,N1] = meshElement(macro1,4);
[Enew3,N3] = meshElement(macro3,4);
cFigure; axisGeom; gpatch(Enew1,N1,'g')
hold on; gpatch(Enew3,N3,'b')
[Elements, Nodes] = mergeMesh(Enew1,N1,Enew3,N3); % combine mesh1 and mesh2  (mesh1+mesh2)
cFigure; axisGeom; gpatch(Elements,Nodes,'r')
% [Enew,N] = meshElement(testMacro,4);
% [Enew,N] = meshElement(testMacro,8);
% [Enew,N] = meshElement(testMacro,3);

%2.2

%% Scratch 

% %Retrieve global coords from local coords
% elem = 3;
% botleft = N(-0.8,-0.8)*V(E{elem},:); %this will work for correctly wound only
% topleft = N(-0.8,0.8)*V(E{elem},:);
% botright = N(0.8,-0.8)*V(E{elem},:);
% topright = N(0.8,0.8)*V(E{elem},:);
% plot(botleft(1),botleft(2),'r.','MarkerSize',markerSize)
% plot(topleft(1),topleft(2),'b.','MarkerSize',markerSize)
% plot(botright(1),botright(2),'g.','MarkerSize',markerSize)
% plot(topright(1),topright(2),'m.','MarkerSize',markerSize)


%% Functions 

%% Attempts at 1

function [E,local_nodes] = innerMesh(mesh_refinement_factor,type)
    %
    if strcmp(type,'biquadratic')
        %Create node coordinates for a uniform 3x3 grid centered at 0
        [Vx,Vy] = ndgrid(-1:1,-1:1);
        V = [Vx(:) Vy(:)];
        
        %Create biquadratic macro element with right-hand node-numbering
        E = [1 3 9 7 2 6 8 4 5];
    elseif strcmp(type,'quadratic triangular')
        %eeeeeeeeeeee
    else
        error('element type not recognized')
    end
    

    
    
end

function [Enew,localNodes] = myInnerMesh(E,V,meshRefinementFactor,type)
    
    if strcmp(type,'biquadratic')
        Enew = [1 5 9
                9 8 1
                5 2 6 
                6 9 5 
                9 6 3 
                3 7 9 
                8 9 7 
                7 4 8];
        
    elseif strcmp(type,'quadratic triangular')
        E = NaN; 
        localNodes = NaN;
    end
end

%% Attempts at 2.2

function [E,N] = meshElement(macro,mesh_refinement_factor)
     %global x,y as fcn of eta, xi = N(eta,xi)*macro


end