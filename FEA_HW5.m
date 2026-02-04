%Ben Sarfati 885573816
%FEA HW5
clear; close all; clc

mrfs = 3;
int = zeros(length(mrfs),1);
for j = 1:length(mrfs)
    %% Parameters 

    %Mesh refinement factor
    mesh_refinement_factor = mrfs(j);

    %Whether or not to flip the clockwise-wound part of the mesh
    correctWinding = true;
    
    %Plot settings 
    fontSize = 30;
    axFontSize = 30;
    lineWidth = 2;
    markerSize = 24;
    textSpacing = 0.1;
    
    %% HW4 stuff
    
    %Create node coordinates for a uniform 3x3 grid centered at 0
    [Vx,Vy] = ndgrid(-1:1,-1:1);
    V = [Vx(:) Vy(:)];
    
    %Create all necessary node coords from the above node coords
    V = [3*V; 5*V([1 3 7 9],:); 7*V([1 2 3 4 6 7 8 9],:); 10*V([2 4 6 8],:)];
    if ~correctWinding
        E = {[1 9 7 5 8 4]... %Left triangle
         [1 9 3 5 6 2]... %Right triangle (wound clockwise)
         [7 9 21 19 8 13 25 12 20]... %Top biquad
         [1 7 19 14 4 12 23 10 17]... %Left biquad
         [3 1 14 16 2 10 22 11 15]... %Bottom biquad
         [9 3 16 21 6 11 24 13 18]}; %Right biquad
    else
        E = {[1 9 7 5 8 4]... %Left triangle
         [9 1 3 5 2 6]... %Right triangle (correctly wound)
         [7 9 21 19 8 13 25 12 20]... %Top biquad
         [1 7 19 14 4 12 23 10 17]... %Left biquad
         [3 1 14 16 2 10 22 11 15]... %Bottom biquad
         [9 3 16 21 6 11 24 13 18]}; %Right biquad
    end
    
    %Mesh elements
    Emeshed = cell(6,1);
    Nmeshed = cell(6,1);
    for i = 1:6
        [Emeshed{i},Nmeshed{i}] = meshElement(V(E{i},:),mesh_refinement_factor);
    end
    
    %Merge meshes
    Emerged = Emeshed{1};
    Nmerged = Nmeshed{1};
    for i = 2:6
        [Emerged,Nmerged] = mergeMesh(Emerged,Nmerged,Emeshed{i},Nmeshed{i});
    end
        
    %% HW5 stuff
    
    %integrate
    int(j) = integrateGeom(ones(size(Nmerged,1),1),Nmerged,Emerged);
end

%% loglog Plot

circleArea = pi*10^2;
errs = abs(int-circleArea)/circleArea;

figure; grid; set(gcf,'color','w');
plot(log10(mrfs),log10(errs),'*r','MarkerSize',markerSize);
xlabel('$\log_{10}(r)$','FontSize',fontSize,'Interpreter','latex')
ylabel('$\log_{10}(e)$','FontSize',fontSize,'Interpreter','latex')
title('Level of Refinement vs. Relative Error of "circle" Area','FontSize',fontSize)
set(gca,'fontSize',fontSize)

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