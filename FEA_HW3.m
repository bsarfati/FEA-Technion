%Ben Sarfati 885573816
%FEA HW3
clear; close all; clc

%% Parameters 

%mesh_refinement_factor
mesh_refinement_factor = 3;

%Plot settings 
fontSize = 30;
axFontSize = 30;
lineWidth = 2;
markerSize = 24;
textSpacing = 0.1;

%% Q1

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


%% Q2

%Create all necessary node coords from the above node coords
V = [3*V; 5*V([1 3 7 9],:); 7*V([1 2 3 4 6 7 8 9],:); 10*V([2 4 6 8],:)];
E = {[1 9 7 5 8 4]... %Left triangle
     [1 9 3 5 6 2]... %Right triangle (wound clockwise)
     [7 9 21 19 8 13 25 12 20]... %Top biquad
     [1 7 19 14 4 12 23 10 17]... %Left biquad
     [3 1 14 16 2 10 22 11 15]... %Bottom biquad
     [9 3 16 21 6 11 24 13 18]}; %Right biquad

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

%Plot everything
subplot(2,1,1);
    hold on;
    title('Global mesh before merge')
    for i = 1:6
        plotElement(Emeshed{i},Nmeshed{i});
    end

subplot(2,1,2);
    hold on
    title('Global mesh after merge')
    plotElement(Emerged,Nmerged);

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
    