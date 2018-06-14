function [sys,x0,str,ts]=Switch_sigma(t,x,u,flag,beta,lambda,DT)
% This S-Function computes the switching signal for switching adaptive control. 
% The input vector u has m inputs that are the estimation errors from a 
% multi-estimator. The parameters are the weighting factor beta 
% (alpha is taken equal to 1), forgetting factor lambda and Dwell-Time DT.


m=4;
x_p=zeros(m,1);

switch flag,
    % Initialization
    case 0,
        sizes = simsizes;
        sizes.NumContStates  = 0;
        sizes.NumDiscStates  = m+2;
        sizes.NumOutputs     = 1;
        sizes.NumInputs      = m;
        sizes.DirFeedthrough = 1;
        sizes.NumSampleTimes = 1;

        sys = simsizes(sizes);
        
        x0  = [zeros(m,1);1;DT];
        str = [];
        ts  = [-1 0]; 
        
     % state update 
    case 2
        % Based on u(i) compute the monintoring signal Ji in a recursive
        % way for all inputs (prediction errors)
        if x(m+2) == DT
            best_Ji = 100000;
            best_i = 1;
            for i=1:m
                error = u(i);
                old_sum = x(m+1);
                new_sum = exp(-lambda)*old_sum + error^2;
                Ji = beta* new_sum + error^2;

                if Ji < best_Ji
                    best_Ji = Ji;
                    best_i = i;
                end

                x(i) = new_sum;            
            end
            x(m+1) = best_i;
            if best_i == 100000
                "ERROR IN SWITCH SIGMA"
            end
        end
        x(m+1) = best_i;
        x(m+2) = x(m+2)-1;
        
        if x(m+2) == 0
            x(m+2) = DT;
        end
        
            
        % Write an algorithm for the Dwell-Time
        % DT shows the number of sampling period of waiting between two
        % switchings and is saved in x(m+2).
        % x(m+1) is the choice of the best predictor.
        
        sys=[x(1:m);x(m+1);x(m+2)];
        
     % output update
    case 3
        disp("Switch")
        disp(x(m+1))
        
        sys=x(m+1); 
    case 9
        sys=[];
 end
    
        
        