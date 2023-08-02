function q = q_vlp(P1, P2, s_P, f_fact, T, Z, theta, D)
q = sqrt((P2^2 - (exp(s_P)*P1^2)) / ((2.685*10^-3*f_fact*Z*T*(exp(s_P)-1))/(sind(theta)*D^5)));
end