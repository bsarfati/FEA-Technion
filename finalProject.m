%Ben Sarfati 885573816
%FEA HW2
clear; close all; clc

%% Parameters

%Mesh refinement factor
mesh_refinement_factor = 1;

%Project problem parameters
a = 2;
b = 5;

%% Mesh project problem 

%Create ordered node coordinates for a single macro element
V = [0 0 
     b 0 
     b a/2
     0 a/2
     b/2 0
     b a/4
     b/2 a/2
     0 a/4
     b/2 a/4];

% % %Full version in case it's not actually sym because stiffness
% % V = [0 0 
% %      b 0 
% %      b a
% %      0 a
% %      b/2 0
% %      b a/2
% %      b/2 a
% %      0 a/2
% %      b/2 a/2];

[Emeshed,Nmeshed] = meshElement(V,mesh_refinement_factor);

%% Vis

figure; hold on;
plotElement(Emeshed,Nmeshed);

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