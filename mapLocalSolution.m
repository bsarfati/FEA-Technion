function [E,N,wVals,gradwVals] = mapLocalSolution(macro,mesh_refinement_factor,toType,aMacro)
% mapLocalSolution - Provides FE solution across an element
%   [E,N] = mapLocalSolution(macro,mesh_refinement_factor,toType,a)
%   IDENTICAL TO meshElement EXCEPT: this function also accepts a, the
%   numerical solution at the macro nodes, and outputs the value of a at
%   local_nodes
% Ben Sarfati 2/2026

%Determine the type of element that requires meshing
    switch size(macro,1)
        case 9
            fromType = 'biquadratic';
        case 6
            fromType = 'quadratic triangular';
        case 3
            fromType = 'linear triangular';
    end

%Retrieve basis functions
[phi,Bhat] = retrieveMapping(fromType);

%Retrieve the mesh at the desired fineness in the reference configuration
[E,local_nodes] = innerMesh(mesh_refinement_factor,fromType,toType);

%Optional: check jacobians
% detJ = @(xi,eta) det(Bhat(xi,eta)*macro);
% for i = 1:length(local_nodes)
%     detJs(i) = detJ(local_nodes(i,1),local_nodes(i,2));
% end
% warning(['smallest Jacobian of transform evaluated at new nodes: ' num2str(min(detJs))]);

%Transform node coordinates from reference configuration into global 
N = zeros(size(local_nodes));
wVals = zeros(size(local_nodes,1),1);
gradwVals = zeros(2,size(local_nodes,1));
for i = 1:size(local_nodes,1)
    N(i,:) = phi(local_nodes(i,1),local_nodes(i,2))*macro;
    wVals(i) = phi(local_nodes(i,1),local_nodes(i,2))*aMacro;

    Je = Bhat(local_nodes(i,1),local_nodes(i,2))*macro;
    B = Je\Bhat(local_nodes(i,1),local_nodes(i,2));
    gradwVals(:,i) = B*aMacro;
end