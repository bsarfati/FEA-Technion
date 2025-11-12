function N = getLagrangePoly(xi)
% getLagrangePoly gets Lagrange polynomials
%   N = getLagrangePoly(xi)
%   returns a function handle containing each of the Lagrange polynomials
%   based on a column vector of nodal position values xi.
%   
% Ben Sarfati 11/2025

N = @(x) arrayfun(@(i) prod(x-xi([1:i-1 i+1:end])',2)/prod(xi(i)-xi([1:i-1 i+1:end])',2), 1:length(xi));

% function N = getLagrangePoly(xi)
% % getLagrangePoly gets Lagrange polynomials
% %   N = getLagrangePoly(xi)
% %   returns a cell array containing function handles of each of the
% %   Lagrange polynomials based on a column vector of nodal position values
% %   xi. Assumes column nodes.
% %   
% % Ben Sarfati 11/2025
% 
% N = cell(length(xi),1);
% for i = 1:length(xi)
%     N{i} = @(x) prod(x-xi([1:i-1 i+1:end])',2)/prod(xi(i)-xi([1:i-1 i+1:end])',2);
% end