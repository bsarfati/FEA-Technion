function [E,local_nodes] = innerMesh(k,type)
% innerMesh - meshes a reference element
%   [E,local_nodes] = innerMesh(k,type)
%   Returns an element list E and node coordinates list local_nodes of a 
%   reference element after subdivision k times on each axis into 2*(k)^2
%   linear triangles.
% Ben Sarfati 1/2026

%Initialize node coordinates and element lists
[Vx,Vy] = ndgrid(linspace(-1,1,k+1),linspace(-1,1,k+1));
local_nodes = [Vx(:) Vy(:)];
E = zeros(2*k^2,3);

%Check requested element type
if strcmp(type,'biquadratic')
    %Create element list
    for i = 1:k
        for j = 1:k
            E(k*(i-1)+j,:) = [(i-1)*(k+1)+j+1 i*(k+1)+j (i-1)*(k+1)+j];
            E(k^2+k*(i-1)+j,:) = [i*(k+1)+j (i-1)*(k+1)+j+1 i*(k+1)+j+1];
        end
    end
elseif strcmp(type,'quadratic triangular')
    %Update local node coordinates list (only have nodes in bottom right)
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
    error('element type not recognized')
end