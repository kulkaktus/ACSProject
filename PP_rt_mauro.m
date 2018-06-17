function [sys,x0,str,ts]=PP_rt(t,x,u,flag,nA,nB,d,P,Hs,Hr,Ts)
% This S-Function computes an RST controller based on the pole placement technique.
% The input vector u has two vector inputs for the coefficients of A and B polynomials 
% and a vector for the desired closed-loop polynomial P. The
% fixed parts of the controllers are specified in Hs and Hr vectors.
% The output of the block is the RST controller as column vectors. 
% The polynomial T is computed to have the same dynamics for regulation and tracking.
% Ts is the sampling period 

nHs=length(Hs)-1;
nHr=length(Hr)-1;

switch flag
    % Initialization
    case 0
        sizes = simsizes;
        sizes.NumContStates  = 0;
        sizes.NumDiscStates  = 0;
        sizes.NumOutputs     = nA+nB+d+2*nHr+2*nHs+1; % thank you for making code working for Hr and Hs arbitrarely sized!
        sizes.NumInputs      = nA+nB;
        sizes.DirFeedthrough = nA+nB+d+nHr+nHs+3;
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
        A = [1;u(1:nA)]';
        B = [zeros(d+1,1);u(nA+1:nA+nB)]';
        
        B = B(d+2:end); % plant model q^-d*B/A

        % orders of polynomials
        nA = length(A) - 1;
        nB = length(B);
        nHr = length(Hr) - 1;
        nHs = length(Hs) - 1;
        nBprime = nB + nHr; % B' = B*Hr
        nAprime = nA + nHs;

        Aprime = conv(A, Hs); % coefficients for polynomial A*Hs
        Bprime = conv(B, Hr);
        
        % Add internal high frequency poles as aux poles to P
        intPoles = roots(A);
        intPoles = intPoles(abs(imag(intPoles)) > 0.001);

        P = conv(P,poly(intPoles));

        alpha = 0.2;
        nP_max = nA + nHs + nB + nHr + d -1;
        n = nP_max - (length(P)-1);
        
        Pf = [1, -alpha];
        for ii = 2:n
            Pf = conv(Pf,[1, -alpha]); 
        end

        P = conv(P, Pf);
        nP = length(P) - 1;

        p = [P, zeros(1,nAprime+nBprime+d-1-nP)]';

        % Construct first part of Sylvester matrix M
        tri_tmp = tril(toeplitz([Aprime,zeros(1,nBprime+d-1)]));
        MA = tri_tmp(:,1:nBprime+d);

        % Construct second part of Sylvester matrix M
        tri_tmp = tril(toeplitz([zeros(1,d+1),Bprime,zeros(1,nAprime-1)]));
        MB = tri_tmp(:,1:nAprime);

        M = [MA,MB];
        
        if rcond(M) > eps % USE A MORE RESTRICTIVE LIMIT?
        
            x = M\p;  % solves A'S' + q^(-d)B'R' = P

            Sprime = x(1:nBprime+d)';
            Rprime = x(nBprime+d+1:end)';
            S = conv(Sprime, Hs);
            R = conv(Rprime, Hr);
        
        else
            % If M is close to singular open the loop by setting R = [0 ... S =
            % [1 0000 ] -> T = 0, u = 0
            R = zeros(1,nAprime+nHr);
            S = [1, zeros(1,nBprime+d+nHs-1)];
            
        end
        
%         % TO AVOID NANs
%         if (any(~isfinite(S)) || any(~isfinite(R)))
%             S = [1, zeros(1,length(S)-1)];
%             R = zeros(1,length(R));
%         end
        
        T = sum(R);

        sys=[R';S';T'];
    case 9
        sys=[];
 end
    
        
        