load('K1.mat');
load('K2.mat');
load('K3.mat');
load('Tormod.mat');
Ts = 0.04;
% S1 = S1';
% S2 = S2';
% R1 = R1';
% R2 = R2';
R1 = [R1,0];
S1 = [S1,0];


A = G3.F;
B = G3.B;
Ts = 0.04;
nA = length(A) - 1;
d = 0;
if length(B) > 1
    for i=2:length(B)
        if B(i) == 0
            d = d + 1;
        end
    end
end
nB = length(B)-d-1;

p1 = -1.835; % rise time 0.7 MAYBE THE POLE VALUES WORK BETTER FOR THE ADAPTIVE ONE
p2 = 0.848; % rise time 0.7
P = [1 p1 p2];

Hs = [1 -1];
Hr = [1 1]; 

beta = 1;
lambda = 0.3;
DT = 25;