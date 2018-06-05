function [R,S]=Poleplace(B,A,Hr,Hs,P)
    %Determine delay d and separate it from B
    d = 0;
    if length(B) > 1
        for i=2:length(B)
            if B(i) == 0
                d = d + 1;
            end
        end
    end
    B = B(1+d:end);
       
    n_a = length(A)-1;
    n_b = length(B)-1;
    n_hr = length(Hr)-1;
    n_hs = length(Hs)-1;
    n_a1 = n_a+n_hs;
    n_b1 = n_b+n_hr;
    n_s1 = n_b+n_hr+d-1;

    A1 = conv(A,Hs);
    B1 = conv(B,Hr);
   
    % M1, M2 are left and right halves of sylvester matrix
    M1 = tril(toeplitz([A1,zeros(1,n_b1+d-1)]));
    M1 = M1(:, 1:n_b1+d);
    M2 = tril(toeplitz([zeros(1,d+1), B1(2:end), zeros(1,n_a1-1)]));
    M2 = M2(:, 1:n_a1);
    
    M = [M1, M2];
    p = ([P,zeros(1, n_a1+n_b1+d-length(P))])';
    x = M\p;
     
    S1 = x(1:n_s1+1)';
    R1 = x(n_s1+2:end)';
    
    S = conv(S1,Hs);
    R = conv(R1,Hr);
end