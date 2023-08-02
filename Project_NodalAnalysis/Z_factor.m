% Beggs and Brill correlation for Z_factor as publischied in J Petrol
% Explor Prod Technol (2016) (Kareem A. L. et al (2016)

function Z = Z_factor(T_pr, P_pr)


A = 1.39.*((T_pr - 0.92).^0.5) - 0.36.*T_pr - 0.1;
E = 9.*(T_pr-1);
B = (0.62 - 0.23.*T_pr).*P_pr + (0.06./(T_pr-0.86) - 0.037).*(P_pr.^2) + (0.32.*P_pr)/10.^E;
C = 0.132 - 0.32.*log10(T_pr);
F = 0.3106 - 0.49.*T_pr + 0.1824.*(T_pr.^2);
D = 10.^F;
Z = A + (1-A)./exp(B) + C.*(P_pr.^D);

end