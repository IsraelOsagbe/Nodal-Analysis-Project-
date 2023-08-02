function P2 = P2_cal_ver(P1, s_P, f_fact, Z, T, Q, theta, D)

 P2 = sqrt((exp(-s_P)*P1^2) - ((2.689e-3*f_fact*Z*T*(Q^2)*(1-exp(-s_P)))/(sind(theta)*D^5)));
 
end