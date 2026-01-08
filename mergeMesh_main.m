close all
clear all
clc


nodes_template = [0 0; 1 0; 1 1; 0 1];
elements_template = [1 2 3;
                    3 4 1];

nodes1 = nodes_template; 
elements1 = elements_template;

nodes2 = nodes_template+[1, 0]; % shift template mesh right
elements2 = elements_template;

nodes3 = nodes_template+[0, 1]; % shift template mesh up
elements3 = elements_template;

nodes4 = nodes_template-[1, 0]; % shift template mesh left
elements4 = elements_template;

nodes5 = nodes_template-[0, 1]; % shift template mesh down
elements5 = elements_template;

[Elements, Nodes] = mergeMesh(elements1,nodes1,elements2,nodes2); % combine mesh1 and mesh2  (mesh1+mesh2)
[Elements, Nodes] = mergeMesh(Elements,Nodes,elements3,nodes3);  % ...(mesh1+mesh2)+mesh3
[Elements, Nodes] = mergeMesh(Elements,Nodes,elements4,nodes4); % ((mesh1+mesh2)+mesh3)+mesh4
[Elements, Nodes] = mergeMesh(Elements,Nodes,elements5,nodes5); 

subplot(2,1,1);
    hold on;
    title('Global mesh before merge')
    plotElement(elements1,nodes1);
    plotElement(elements2,nodes2);
    plotElement(elements3,nodes3);
    plotElement(elements4,nodes5);
    plotElement(elements5,nodes4);

subplot(2,1,2);
    hold on
    title('Global mesh after merge')
    plotElement(Elements,Nodes);

function [Elements, Nodes] = mergeMesh(E1,N1,E2,N2)
    % E1 - element connectivity table of mesh 1 (current, doesn't change)
    % N1 - nodal coordinates table of mesh 1 (current, doesn't change)      
    % E2 - element connectivity table of mesh 2 (additional)
    % N2 - nodal coordinates table of mesh 2 (additional)      

    eps = 1e-5; % define tolerance distinguishing duplicated nodes (robust to rounding errors)

    n_current = size(N1,1); % number of nodes in current mesh
    E2_shifted = E2+n_current; % offset nodal ID's of mesh 2 by n_current (assume no dup nodes)
    new_nodes = N2;
    for ii=1:size(N2)
        temp_node_pos = N2(ii,:);
        % determine if temp_node_pos is a dup node (search can be
        % restricted to boundary nodes for efficiency)
        d=((N1(:,1)-temp_node_pos(1)).^2+(N1(:,2)-temp_node_pos(2)).^2).^.5; % distances to all current nodes
        [min_distance, min_ind] = min(d); % find nearest neighbor 
        if min_distance<=eps % check coindidence
            % node ii from mesh 2 already appears in N1
            above_ii_mask = E2>ii;
            equal_ii_mask = E2==ii;
            E2_shifted(above_ii_mask) = E2_shifted(above_ii_mask)-1; % shift node id's 1 down (only from ii+1 and above)
            E2_shifted(equal_ii_mask) = min_ind; % change all ocurances of ii in E2 to the id as it appears in the current mesh
            new_nodes(ii,1) = NaN; % mark row for deletion
        end
    end
    new_nodes(isnan(new_nodes(:,1)),:) = [];
    Elements = [E1; E2_shifted];
    Nodes = [N1; new_nodes];
end



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
    



