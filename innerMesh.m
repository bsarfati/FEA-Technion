function [E,local_nodes] = innerMesh(k,fromType,toType)
% innerMesh - meshes a reference element
%   [E,local_nodes] = innerMesh(k,type)
%   Returns an element list E and node coordinates list local_nodes of a 
%   reference element after subdivision k times on each axis into 2*(k)^2
%   triangular elements of type toType.
% Ben Sarfati 2/2026

%Check requested macro element type
if strcmp(fromType,'biquadratic')
    %Check requested output element type
    if strcmp(toType,'linear triangular')
        %Create node coordinates list
        [Vx,Vy] = ndgrid(linspace(-1,1,k+1),linspace(-1,1,k+1));
        local_nodes = [Vx(:) Vy(:)];
    
        %Create element list
        E = zeros(2*k^2,3);
        for i = 1:k
            for j = 1:k
                E(k*(i-1)+j,:) = [(i-1)*(k+1)+j+1 i*(k+1)+j (i-1)*(k+1)+j];
                E(k^2+k*(i-1)+j,:) = [i*(k+1)+j (i-1)*(k+1)+j+1 i*(k+1)+j+1];
            end
        end
    else
        if mod(k,2)
            error('Selected number of mesh refinements is odd; quadratic triangles cannot be created.')
        end

        %Create node coordinates list
        [Vx,Vy] = ndgrid(linspace(-1,1,2*k+1),linspace(-1,1,2*k+1));
        local_nodes = [Vx(:) Vy(:)];
    
        %Create element list
        E = zeros(2*k^2,6);
        for i = 1:2:2*k-1
            for j = 1:2:2*k-1
                E(k*((i+1)/2-1)+(j+1)/2,:) = [(i-1)*(2*k+1)+j+2 (i+1)*(2*k+1)+j (i-1)*(2*k+1)+j i*(2*k+1)+j+1 i*(2*k+1)+j (i-1)*(2*k+1)+j+1];
                E(k^2+k*((i+1)/2-1)+(j+1)/2,:) = [(i+1)*(2*k+1)+j (i-1)*(2*k+1)+j+2 (i+1)*(2*k+1)+j+2 i*(2*k+1)+j+1 i*(2*k+1)+j+2 (i+1)*(2*k+1)+j+1];
            end
        end
    end
elseif strcmp(fromType,'quadratic triangular') | strcmp(fromType,'linear triangular')
    %Check requested output element type
    if strcmp(toType,'linear triangular')
        %Create local node coordinates list (only have nodes in bottom left)
        [Vx,Vy] = ndgrid(linspace(0,1,k+1),linspace(0,1,k+1));
        local_nodes = [Vx(:) Vy(:)];
        local_nodes = local_nodes(logical(reshape(flipud(tril(ones(k+1))),[],1)),:);
        
        %Create element list
        % E = zeros(2*k^2,3);
        E = [];
        for i = 1:k
            for j = 1:k+1-i
                %Complicated triangle version figured out on paper
                E(end+1,:) = [1/2*(i-1)*(2-i)+(i-1)*(k+1)+j+1 1/2*i*(1-i)+i*(k+1)+j 1/2*(i-1)*(2-i)+(i-1)*(k+1)+j];
                % E(end+1,:) = [1/2*(i-1)*(2-i)+(i-1)*(k+1)+j 1/2*(i-1)*(2-i)+(i-1)*(k+1)+j+1 1/2*i*(1-i)+i*(k+1)+j];
                if j ~= 1
                    E(end+1,:) = [1/2*i*(1-i)+i*(k+1)+j-1 1/2*(i-1)*(2-i)+(i-1)*(k+1)+j 1/2*i*(1-i)+i*(k+1)+j];
                    % E(k*(k+1)/2+,:) = [1/2*(i-1)*(2-i)+i*(k+1)+j+1 1/2*(i-1)*(2-i)+i*(k+1)+j 1/2*(i-1)*(2-i)+(i-1)*(k+1)+j];
                end
            end
        end
    else 
        %For quadratic tri to quadratic tri the math is even harder. We can
        %just chop the biquadratic output in 2.
        if mod(k,2)
            error('Selected number of mesh refinements is odd; quadratic triangles cannot be created.')
        end

        %Create node coordinates list
        [Vx,Vy] = ndgrid(linspace(0,1,2*k+1),linspace(0,1,2*k+1));
        local_nodesRaw = [Vx(:) Vy(:)];
    
        %Create element list
        Eraw = zeros(2*k^2,6);
        for i = 1:2:2*k-1
            for j = 1:2:2*k-1
                Eraw(k*((i+1)/2-1)+(j+1)/2,:) = [(i-1)*(2*k+1)+j+2 (i+1)*(2*k+1)+j (i-1)*(2*k+1)+j i*(2*k+1)+j+1 i*(2*k+1)+j (i-1)*(2*k+1)+j+1];
                Eraw(k^2+k*((i+1)/2-1)+(j+1)/2,:) = [(i+1)*(2*k+1)+j (i-1)*(2*k+1)+j+2 (i+1)*(2*k+1)+j+2 i*(2*k+1)+j+1 i*(2*k+1)+j+2 (i+1)*(2*k+1)+j+1];
            end
        end

        %Find all nodes that belong to the lower left triangle
        logicalKeepNodes = logical(reshape(flipud(tril(ones(2*k+1))),[],1));

        %Find new numberings of nodes in E
        [~,newElemNumbering] = ismember(Eraw,find(logicalKeepNodes));
        
        %Remove elements that contained nodes not belonging to subset
        E = newElemNumbering(all(newElemNumbering~=0,2),:);

        %Update node list accordingly
        local_nodes = local_nodesRaw(logicalKeepNodes,:);

    end
else
    error('element type not recognized')
end