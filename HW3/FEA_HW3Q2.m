%% FEA HW3 Q2
%Ben Sarfati 941180069

clear all; close all; clc;

%%

E = 1E10; A = 0.0001; L = 1; sigmay = 276E6; k0 = 1E50; P = 1E4;
% E = 1E10; A = 0.0001; L = 1; sigmay = 276E6; k0 = 1E50; P = 39837;

Khat = E*A/L*[1 -1; -1 1];
alpha = [pi/3 -pi/3 0];
T = @(e) [cos(alpha(e)) sin(alpha(e)) 0 0; 0 0 cos(alpha(e)) sin(alpha(e))];

Kg_expanded = zeros(12);
for e = 1:3
    K = T(e)'*Khat*T(e);
    Kg_expanded(2*e-1:2*e+2,2*e-1:2*e+2) = Kg_expanded(2*e-1:2*e+2,2*e-1:2*e+2) + K;
end
Kg = Kg_expanded(1:6,1:6)+Kg_expanded(7:12,1:6)+Kg_expanded(1:6,7:12)+Kg_expanded(7:12,7:12);

%%Boundary conditions and accounting for F
Kg(1,1) = Kg(1,1) + k0;
Kg(2,2) = Kg(2,2) + k0;
Kg(5,5) = Kg(5,5) + k0;
Kg(6,6) = Kg(6,6) + k0;
Kg(4,:) = Kg(4,:) +1/sqrt(3)*Kg(3,:);
row3 = Kg(3,:);
Kg(3,:) = [0 0 sqrt(3) -1 0 0];
Fg = [0 0 0 P 0 0]';

a = Kg\Fg;

strain = -diff(T(1)*a(3:6))/L;
stress = E*strain;
stress/sigmay