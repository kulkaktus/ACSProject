function [sys,x0,str,ts]=Switch_sigma(t,x,u,flag,beta,lambda,DT)
% This S-Function computes the switching signal for switching adaptive control.
% The input vector u has m inputs that are the estimation errors from a
% multi-estimator. The parameters are the weighting factor beta
% (alpha is taken equal to 1), forgetting factor lambda and Dwell-Time DT.


m=4;
x_p=zeros(m,1);
alpha = 1;

switch flag
    % Initialization
    case 0
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
        Ji = zeros(1,m);
        prev = zeros(1,m);
        for i=1:m
            Ji(i) = (alpha+beta)*u(i)^2 + x(i);
            prev(i) = exp(-lambda)*(Ji(i) - alpha*u(i)^2);
        end
        
        % Write an algorithm for the Dwell-Time
        % DT shows the number of sampling period of waiting between two
        % switchings and is saved in x(m+2).
        % x(m+1) is the choice of the best predictor.
        [~, sigma] = min(Ji);
        
        % Dwell-time (simple waiting)
        
%         if (sigma ~= x(m+1))
%             if (x(m+2) > DT)
%                 x(m+2) = 0;
%                 x(m+1) = sigma;
%             else
%                 x(m+2) = x(m+2) + 1;
%             end
%             
%         end 

        % allowed to switch after DT sampling periods
        % allows faster switching when sigma stayed constant over more than
        % 1 DT
        if ((sigma ~= x(m+1)) && (x(m+2) > DT))
            x(m+2) = 0;
            x(m+1) = sigma;
        else
            x(m+2) = x(m+2) + 1;
        end
        % save prev, sigma, passed time since last switching, minimal cost
        sys=[prev';x(m+1);x(m+2)];
        
        % output update
    case 3
        sys=x(m+1);
    case 9
        sys=[];
end


