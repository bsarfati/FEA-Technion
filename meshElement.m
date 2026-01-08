function [E,N] = meshElement(macro,mesh_refinement_factor)
% innerMesh - meshes a reference element
%   [E,N] = meshElement(macro,mesh_refinement_factor)
%   Meshes an element with the node coordinates list given by macro into 
%   linear triangular elements, subdividing into mesh_refinement_factor for
%   each side
% Ben Sarfati 1/2026

%Determine the type of element that requires meshing
if size(macro,1) == 9
    type = 'biquadratic';

    %Define isoparametric mapping
    phi = @(xi,eta) (1/4)*[
    xi*(xi-1)*eta*(eta-1)
    xi*(xi+1)*eta*(eta-1)
    xi*(xi+1)*eta*(eta+1)
    xi*(xi-1)*eta*(eta+1)
    2*(1-xi^2)*eta*(eta-1)
    2*xi*(xi+1)*(1-eta^2)
    2*(1-xi^2)*eta*(eta+1)
    2*xi*(xi-1)*(1-eta^2)
    4*(1-xi^2)*(1-eta^2)]';
else
    type = 'quadratic triangular';

    %Define isoparametric mapping
    phi = @(xi,eta) [
    xi*(2*xi-1)
    eta*(2*eta-1)
    (1-xi-eta)*(1-2*xi-2*eta)
    4*xi*eta
    4*eta*(1-xi-eta)
    4*xi*(1-xi-eta)]';
end

%Retrieve the mesh at the desired fineness in the reference configuration
[E,local_nodes] = innerMesh(mesh_refinement_factor,type);

%Transform node coordinates from reference configuration into global 
N = zeros(size(local_nodes));
for i = 1:size(local_nodes,1)
    N(i,:) = phi(local_nodes(i,1),local_nodes(i,2))*macro;
end