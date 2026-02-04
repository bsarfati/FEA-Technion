function int = integrateGeom(u,N,E)
% integrateGeom - performs numerical integration over a meshed surface
%   int = integrateGeom(u,N,E)
%   returns the result of numerical integration of a function given by
%   nodal values u at nodal coordinates N meshed as linear triangular
%   elements listed in E.
% Ben Sarfati 1/2026

%Calculate individual contribution of each element’s integral to the total
int = 0;
for i = 1:size(E,1)
	%Retrieve value of function at element’s nodes
	uNodal = u(E(i,:));

    %Retrieve global positions of local nodes
    macro = N(E(i,:),:);

	%Perform gauss quadrature in reference configuration (1 point exact;
	% reduces to area of element * average of nodal values; area = 0.5)
	localInt = 1/6*(sum(uNodal));

    %Retrieve gradient of isoparametric mapping with respect to reference 
    % configuration for linear triangular elements
    [~,Bhat] = retrieveMapping('linear triangular');

	%Compute determinant of jacobian of current element
	J = Bhat*macro;
    detJ = det(J);

	%Add element’s contribution to total integral value
	int = int+detJ*localInt;
end
