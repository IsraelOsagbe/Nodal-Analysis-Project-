function Z = Z_factor2(z, rho)

A1 = 0.3265;
A2 = -1.0700;
A3 = -0.5339;
A4 = 0.01569;
A5 = -0.05165;
A6 = 0.5475;
A7 = -0.7361;
A8 = 0.1844;
A9 = 0.1056;
A10 = 0.6134;
A11 = 0.7210;
z = 1;

rho = (0.27*P_pr)/(z*T_pr);
Z = 1 + rho*(A1 + A2/T_pr + A3/T_pr^3 + A4/T_pr^4 + A5/T_pr^5) + rho^2*(A6 + A7/T_pr + A8/T_pr^2) - A9*rho^5*(A7/T_pr + A8/T_pr^2) + A10*(1 + A11*rho^2)*(rho^2/T_pr^3)*exp(-A11*rho^2);

end