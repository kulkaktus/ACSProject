Ts= 0.04;
M_p = 1.05;
T_r = 0.6;
zeta = 1/sqrt(1+pi^2/(log(M_p-1)^2));


omega_n = (2.16*zeta+0.6)/T_r
p1 = -2*exp(-zeta*omega_n*Ts)*cos(omega_n*Ts*sqrt(1-zeta^2))
p2 = exp(-2*zeta*omega_n*Ts)