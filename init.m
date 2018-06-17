load('K1.mat');
load('K2.mat');
load('K3_with_four_poles.mat');
load('Tormod.mat');

R3 = [R3,0];
S3 = [S3,0];

Ts = 0.04;
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

p1 = -1.8555; % rise time 0.7 MAYBE THE POLE VALUES WORK BETTER FOR THE ADAPTIVE ONE
p2 = 0.8656; % rise time 0.7
P = [1 p1 p2];

Hs = [1 -1];
Hr = [1 1]; 

beta = 40; 
lambda = 0.95; %0.2-0.6
DT = 25;