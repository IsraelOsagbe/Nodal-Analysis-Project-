function P = P2_cal(P, s_P, f_fact, Z, T, Q, theta, D)

 P = sqrt((exp(s_P)*P^2) + ((2.689*10^-3*f_fact*Z*T*(Q^2)*(exp(s_P) - 1))/(sind(theta)*D^5)));
 
end