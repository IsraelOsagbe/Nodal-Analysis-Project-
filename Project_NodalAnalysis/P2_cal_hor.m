function P2 = P2_cal_hor(P1, Y_g, f_fact, Z, T, Q, D, L)

 P2 = sqrt((P1^2) + ((1.007e-4*Y_g*f_fact*Z*T*(Q^2)*L)/(D^5)));
 
end