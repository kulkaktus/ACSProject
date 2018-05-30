function [sys,x0,str,ts]=LPV_Torsional(t,x,u,flag)
% This S-Function computes the output signal of the torsional system at 
% different operating point defined by a scheduling parameter theta. 
% The input vector u includes the input signal of the plant model and the
% scheduling parameter alpha. 
% The output of the block is the output signal. The vector x includes the 
% state of the system (passed inputs and outputs).


% The model parameters for theta=-1
B1=[0  0 -0.006243175204197 0.011753895243426 -0.005224306511396 0.006818043079367
];
A1=[1.0000   -4.244157895679979    8.478830503146163  -10.376345593429134    8.217248271243790 ...
    -3.983451401674943    0.907953679231588];

% The model parameters for theta=0
B2=[0         0    0.017887692424122   -0.041030304411798    0.047476716374229   -0.007411346884422];
A2=[1.0000   -3.960284879882827    7.510091428650813   -8.877259828787567    6.906141956331947 ...
    -3.387645223289709    0.808986771771818];

% The model parameters for theta=1
B3=[0         0    0.008025961669808   -0.011836509528379    0.020387942539163    0.007898229769799];
A3=[1.0000   -3.734863259924208    6.958900190116942   -8.340738205000671    6.680557173866394 ...
    -3.438590193475025    0.874985026452499];

Ts=0.04;

nA=length(A1)-1;
nB=length(B1)-1;

n=nA+nB;

switch flag,
    % Initialization
    case 0,
        sizes = simsizes;
        sizes.NumContStates  = 0;
        sizes.NumDiscStates  = n+1;
        sizes.NumOutputs     = 1;
        sizes.NumInputs      = 2;
        sizes.DirFeedthrough = 0;
        sizes.NumSampleTimes = 1;

        sys = simsizes(sizes);

        x0  = zeros(sizes.NumDiscStates,1);
        str = [];
        ts  = [Ts 0]; 
        
     % state update 
    case 2
        u_k=u(1);
        alpha=u(2);
        
        if u(2) < -1, alpha=-1;end
        if u(2) > 1, alpha=1;end
        
        if alpha <= 0,
            B=B2+alpha*(B2-B1);
            A=A2+alpha*(A2-A1);
        else
            B=B2+alpha*(B3-B2);
            A=A2+alpha*(A3-A2);
        end
        
        y_k=[A(2:end) B(2:end)]*x(1:n);

        sys=[-y_k;x(1:nA-1);u(1);x(nA+1:n-1);alpha];
               
     % output update
    case 3  
        alpha=x(end);      
        if alpha <= 0
            B=B2+alpha*(B2-B1);
            A=A2+alpha*(A2-A1);
        else
            B=B2+alpha*(B3-B2);
            A=A2+alpha*(A3-A2);
        end
        
        sys=[A(2:end) B(2:end)]*x(1:n);
         
    case 9
        sys=[];
 end
    