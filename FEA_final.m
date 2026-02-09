%Ben Sarfati 885573816
%FEA HW2
clear; close all; clc

%% Parameters

%Choose geometry
geometries = {'benchmark','project','simple'};
geometry = geometries{1};

%Mesh refinement factor
mesh_refinement_factor = 3;

%Benchmark problem parameters
R = 10; %Benchmark problem radius
% % % % meshBoundarRes = 100; %Number of points on boundary of mesh

%Project problem parameters
a = 2;
b = 5;

%% Mesh benchmark problem (using HW3 geometry)

% % % % thetaVals = linspace(0,2*pi,meshBoundaryRes)

%Create node coordinates for a uniform 3x3 grid centered at 0
[Vx,Vy] = ndgrid(-1:1,-1:1);
V1 = [Vx(:) Vy(:)];

%Create all necessary node coords from the above node coords
V1 = [3*V1; 5*V1([1 3 7 9],:); sqrt(R^2/2)*V1([1 2 3 4 6 7 8 9],:); R*V1([2 4 6 8],:)];
E1 = {[1 9 7 5 8 4]... %Left triangle
 [9 1 3 5 2 6]... %Right triangle (correctly wound)
 [7 9 21 19 8 13 25 12 20]... %Top biquad
 [1 7 19 14 4 12 23 10 17]... %Left biquad
 [3 1 14 16 2 10 22 11 15]... %Bottom biquad
 [9 3 16 21 6 11 24 13 18]}; %Right biquad

%Mesh elements
EmeshedBM = cell(6,1);
NmeshedBM = cell(6,1);
for i = 1:6
    [EmeshedBM{i},NmeshedBM{i}] = meshElement(V1(E1{i},:),mesh_refinement_factor);
end

%Merge meshes
EmergedBM = EmeshedBM{1};
NmergedBM = NmeshedBM{1};
for i = 2:6
    [EmergedBM,NmergedBM] = mergeMesh(EmergedBM,NmergedBM,EmeshedBM{i},NmeshedBM{i});
end
        
%% Mesh project problem 

%Create ordered node coordinates for a single macro element
V2 = [0 0 
     b 0 
     b a/2
     0 a/2
     b/2 0
     b a/4
     b/2 a/2
     0 a/4
     b/2 a/4];

[EmeshedP,NmeshedP] = meshElement(V2,mesh_refinement_factor);

%% Vis
figure; hold on;
plotElement(EmergedBM,NmergedBM);
figure; hold on;
plotElement(EmeshedP,NmeshedP);