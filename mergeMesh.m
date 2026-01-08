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