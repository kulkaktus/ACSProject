function [sys,x0,str,ts]=PP_rt(t,x,u,flag,nA,nB,d,P,Hs,Hr,Ts)
% This S-Function computes an RST controller based on the pole placement technique.
% The input vector u has two vector inputs for the coefficients of A and B polynomials 
% and a vector for the desired closed-loop polynomial P. The
% fixed parts of the controllers are specified in Hs and Hr vectors.
% The output of the block is the RST controller as column vectors. 
% The polynomial T is computed to have the same dynamics for regulation and tracking.
% Ts is the sampling period 

n_Hs=length(Hs)-1;
n_Hr=length(Hr)-1;

switch flag,
    % Initialization
    case 0,
        sizes = simsizes;
        sizes.NumContStates  = 0;
        sizes.NumDiscStates  = 0;
        sizes.NumOutputs     = nA+nB+d+n_Hr+n_Hs+3;
        sizes.NumInputs      = nA+nB;
        sizes.DirFeedthrough = nA+nB+d+n_Hr+n_Hs+3;
        sizes.NumSampleTimes = 1;

        sys = simsizes(sizes);

        x0  = [];
        str = [];
        ts  = [Ts 0]; 
        
     % state update 
    case 2       
        sys=[];
        
     % output update
    case 3
        A=[1;u(1:nA)]';
        B=[zeros(d+1,1);u(nA+1:nB+nA)]';
        
        % The students should write the code for computing an RST
        % controller
        [R,S] = Poleplace(B,A,Hr,Hs,P);
        T = sum(R); %Same tracking and regulation dynamics
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% Should we use Q-parametrization? Too slow here probably %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sys=[R;S;T];
    case 9
        sys=[];
 end
    
        
        