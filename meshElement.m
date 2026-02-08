function [E,N] = meshElement(macro,mesh_refinement_factor)
% innerMesh - meshes a reference element
%   [E,N] = meshElement(macro,mesh_refinement_factor)
%   Meshes an element with the node coordinates list given by macro into 
%   linear triangular elements, subdividing into mesh_refinement_factor for
%   each side
% Ben Sarfati 2/2026

%Determine the type of element that requires meshing
    switch size(macro,1)
        case 9
            type = 'biquadratic';
        case 6
            type = 'quadratic triangular';
        otherwise
            error('element type not implemented');
    end

%Retrieve basis functions
phi = retrieveMapping(type);

%Retrieve the mesh at the desired fineness in the reference configuration
[E,local_nodes] = innerMesh(mesh_refinement_factor,type);

%Transform node coordinates from reference configuration into global 
N = zeros(size(local_nodes));
for i = 1:size(local_nodes,1)
    N(i,:) = phi(local_nodes(i,1),local_nodes(i,2))*macro;
end