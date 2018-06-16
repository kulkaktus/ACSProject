function [sys,x0,str,ts,simStateCompliance]=aap(t,x,u,flag,nA,nB,d,variance_noise)
% This function is a parametric adaptation algorithm. It takes the input
% and output of a system and computes the model parameters. The adaptation 
% gain is initialized to F0 and the parameters are initialized to the
% parameters of the last model chosen by the switching algorithm. The
% function has three inputs: u(1) is the plant input, u(2) the plant output
% and u(3) is the switching signal belongs to {1,2,3,4}. The output of the
% block is the predicted output yhat and the vector of estimated
% parameters theta. The parameter vector theta is initialized by the best 
% model among the fixed one after each transition of the switching signal u(3). 

n=nA+nB;
Nstate=n*(n+2)+d+1;% The last state concernes the value of u(3) in the previous sampling time.
deadzone=2*sqrt(variance_noise); % this value can be chosen to freeze adaptation if the adaptation error is too small.
F0=eye(n);
%tr_F0 = 20*n;
lambda2 = 1;

% The parameters of G1
B1=[0         0   -0.0062    0.0118   -0.0052    0.0068]';
A1=[1.0000   -4.2442    8.4788  -10.3763    8.2172   -3.9835    0.9080]';

% The parameters of G2
B2=[0         0    0.0179   -0.0410    0.0475   -0.0074]';
A2=[1.0000   -3.9603    7.5101   -8.8773    6.9061   -3.3876    0.8090]';

% The parameters of G3
B3=[0         0    0.0080   -0.0118    0.0204    0.0079]';
A3=[1.0000   -3.7349    6.9589   -8.3407    6.6806   -3.4386    0.8750]';



switch flag
    case 0
        
        sizes = simsizes;
        sizes.NumContStates  = 0;
        sizes.NumDiscStates  = Nstate;
        sizes.NumOutputs     = n+1;
        sizes.NumInputs      = 3;
        sizes.DirFeedthrough = 0;
        sizes.NumSampleTimes = 1;

        sys = simsizes(sizes);

        
        x0  = zeros(Nstate,1);
        str = [];
        ts  = [-1 0]; 
        simStateCompliance = 'DefaultSimState';
        
      % state update  
    case 2
        % page 147 script
        % theta_k = [a1,..., an, b1,..., bn]
        % phi_k = [-y(k),..., -y(k-na+1), u(k-d),..., u(k-d-nb+1)]
        
        theta_k = x(1:n);
        phi_k = [x(n+1:n+nA);x(n+nA+d+1:2*n+d)];        
        F_k = reshape(x(2*n+d+1:end-1),n,n);
        
        % Compute the a priori prediction error
        epsilon = u(2) - theta_k'*phi_k;
        
        % update forgetting factor
        alpha = 0.05;
        lambda1 = 1 - alpha*(((epsilon/10))^2/(1+phi_k'*F_k*phi_k));
        %lambda1 = 0.99;
        
        % constant trace
        %alpha = 0.99;
        %lambda1 = max(trace(F_k - (F_k*phi_k*phi_k'*F_k)./(alpha + phi_k'*F_k*phi_k))/tr_F0,0.1);
        %lambda2 = lambda1/alpha;
        
        % Compute F(t+1) using matrix inversion lemma
        F_p = (F_k - (F_k*phi_k*phi_k'*F_k)./(lambda1/lambda2 + phi_k'*F_k*phi_k))./lambda1;
        

        % Dead zone
        if abs(epsilon) >= deadzone
            % Compute F(t+1) using matrix inversion lemma
            %F_p = (F_k - (F_k*phi_k*phi_k'*F_k)./(lambda1/lambda2 + phi_k'*F_k*phi_k))./lambda1;
            % Compute the new estimate for theta
            theta_p = theta_k + F_p*phi_k*epsilon;
            
        else % don't update theta is the adaptation error is small enough (usefull if noise is present)
            theta_p = theta_k;
            %F_p = F_k;
        end
        
        % Update the augmented observation vector
        phi_p = [-u(2);x(n+1:n+nA-1);u(1);x(n+nA+1:2*n+d-1)];
        
        % Initialize the parameters with the best model in the transition 
        % of the switching signal (its past value is saved in x(end))
        
        if u(3)~= x(end)
            % reinitialize the parameters acording to u(3)
            % Which parts do we keep in B1?????
            %F_p = 100*F0;
            switch u(3)
                case 1
                    theta_p = [A1(2:end);B1(d+2:end)]; 
                    F_p = 100*F0;
                    %F_p = tr_F0*F0;
                case 2
                    theta_p = [A2(2:end);B2(d+2:end)];
                    F_p = 100*F0;
                    %F_p = tr_F0*F0;
                case 3
                    theta_p = [A3(2:end);B3(d+2:end)];
                    F_p = 100*F0;
                    %F_p = tr_F0*F0;
                case 4
                    %theta_p = zeros(size(theta_p));
                    %F_p = 100*F0;
            end

        end
        sys = [theta_p;phi_p;reshape(F_p,n*n,1);u(3)];
        
      % output update
    case 3
        % Using strategy 2: sample the output and inputs, compute the
        % controller signal based on the previous estimation of theta (this function outputs theta,
        % which is used by the pole place function to compute the controller parameters and the control signal), then
        % update the estimates of theta (done here in case 2) and the
        % controller parameters based on the new measurements
        theta_p = x(1:n);
        phi_k = [x(n+1:n+nA);x(n+nA+d+1:2*n+d)];
        % Compute yhat and theta_k
        theta_k = theta_p;
        yhat = theta_p'*phi_k;
        sys = [yhat;theta_k];

    case 9
        sys = [];
 end
    
        
        