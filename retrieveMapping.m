function [phi,Bhat] = retrieveMapping(elementType)
% retrieveMapping - retrieves basis functions and derivatives
%   [phi,Bhat] = retrieveMapping(elementType)
%   retrieves isoparametric mapping phi and its gradient Bhat with respect
%   to reference coordinates for various element types
% Ben Sarfati 2/2026

%Determine the type of element that requires meshing
if strcmp(elementType,'biquadratic')
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

    Bhat = NaN; %Needs to be implemented
elseif strcmp(elementType,'quadratic triangular')
    phi = @(xi,eta) [
    xi*(2*xi-1)
    eta*(2*eta-1)
    (1-xi-eta)*(1-2*xi-2*eta)
    4*xi*eta
    4*eta*(1-xi-eta)
    4*xi*(1-xi-eta)]';

    Bhat = NaN; %Needs to be implemented
elseif strcmp(elementType,'linear triangular')
    phi = @(xi,eta) [
    xi
    eta
    1 - xi - eta]';

    Bhat = @(xi,eta) [eye(2) -ones(2,1)];
    % % Bhat = [eye(2) -ones(2,1)]; %Functional for now (FOR INTEGRATEGEOM); 
    % % % may need to become a function handle in the future
end