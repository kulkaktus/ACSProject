load('TorMod.mat');
n_models = 3;
Gs = [G1; G2; G3];
As = [G1.f; G2.f; G3.f];
Bs = [G1.b; G2.b; G3.b];

p1 = -1.7689; 
p2 = 0.7939;
P = [1 p1 p2];

Hs = [1 -1];
Hr = [1 1]; 
a = 0.2;
nq = 9;


["Model", "Margin", "RiseTime", "a", "nq"] 
Candidates = zeros(70, 5); % model, number of candidates, [margin, risetime, a, nq]
index = 1;  
for model=1 %For model 1, a=0.2, nq=9 gives 0.37 MM
    
    for a = 0.27:0.005:0.29 %for 1, tried 0.1:0.02:0.5, none of them are feasible
        for nq = 8:12
            G = Gs(model);
            B = Bs(model,:);
            A = As(model,:);
            poles_aux = [a,a,a,a,a,a,a,a,a,a];
            coefs = poly(poles_aux);
            P_new = conv(P, coefs);

            [R, S] = Poleplace(B, A, Hr, Hs, P_new);
            P_end = conv(A,S) + conv(B,R);

            B = B(2:end); %Separate delay and B
            q_delay = [0 1];

            M_m = 0.4;
            U_max = 56.2; %35 dB  = 56.2

            R = [R, zeros(1, nq)];
            S = [S, zeros(1, nq)];

            R_new = @(Q) R + conv(A, conv(Hr, conv(Hs, Q)));
            S_new = @(Q) (S - conv(q_delay, conv(B, conv(Hs, conv(Hr, Q)))));

            K_new = @(Q) tf(R_new(Q),S_new(Q), Ts,'variable','z^-1');
            output_sensitivity = @(Q) feedback(1,K_new(Q)*G);
            Mod_marg = @(Q) norm(output_sensitivity(Q), Inf)^(-1);

            %Set inequality (c) and equality (ceq) constraints
            c = @(Q) [norm(M_m*output_sensitivity(Q), Inf) - 1;
                 norm(tf(conv(A, R_new(Q)), P_end, Ts, 'variable', 'z^-1'), Inf) - U_max];
            ceq = @(Q) [];
            Nonlincon = @(Q)deal(c(Q), ceq(Q));


            Q_opt = fmincon(Mod_marg, zeros(1, nq),[],[],[],[],[-Inf,-Inf],[Inf,Inf], Nonlincon);
            
            R_final = R_new(Q_opt);
            S_final = S_new(Q_opt);
            T = sum(R_final);

            K_final = tf(R_final,S_final, Ts,'variable','z^-1');
            Sens_out_new = feedback(1,K_final*G);
            MM_new = Mod_marg(Q_opt);
            risetime = 5
            %if MM_new > 0.4
                CL = tf(conv(T,G.b), P_end, Ts,'variable','z^-1'); 
                CL_info = stepinfo(CL);
                risetime = CL_info.RiseTime;
            %end
            Candidates(index,:) = [model, MM_new, risetime, a, nq];
            index = index + 1
            
            
        end
    end
end
"Model candidates:"
Candidates
%% Check constraints
CL = tf(conv(T,G.b), P_end, Ts,'variable','z^-1'); 
U=tf(conv(A,R_final), P_end, Ts, 'variable', 'z^-1');

input= tf(conv(T,A), conv(A,S_final) + conv(conv(q_delay, B), R_final), Ts,'variable','z^-1'); figure(1);  
step(input); 
stepinfo(CL) 
%figure(3);
%bodemag(U);


%% Plots
hold on
figure(1)
subplot(2,2,1)
CL_old = tf(conv(T_old,G.b), P_end_old, Ts,'variable','z^-1'); 
CL = tf(conv(T,G.b), P_end, Ts,'variable','z^-1'); 
title("Tracking step response, original and improved controller")
xlabel( "Time [s]")
legend("Original", "Improved")
step(CL_old, CL)


subplot(2,2,2)
B = G.B;
B = B(2:end);


title("Control signal, original and improved controller")
xlabel("Time [s]")
ylabel("u")
legend("Original", "Improved")
step(input_old, input)


subplot(2,2,3)
Sens_out = feedback(1,K*G);
K_final = tf(R_final,S_final, Ts,'variable','z^-1');
Sens_out_new = feedback(1,K_final*G);
MM_old = norm(Sens_out, Inf)^-1
MM_new = (norm(Sens_out_new, Inf))^(-1)

title("Output sensitivity function, original and improved controller")
legend("Original","Improved")
bodemag(Sens_out, Sens_out_new)


subplot(2,2,4)
U=tf(conv(A,R_old), P_end_old, Ts, 'variable', 'z^-1');
U_new=tf(conv(A,R_final), P_end, Ts, 'variable', 'z^-1');
title("Input sensitivity function, original and improved controller")
legend("Original","Improved")
bodemag(U, U_new)