function [P2, s_P, f_fact, N_Re, m_g, rho_g, B_g, Z_avg] = ultim_func(P1, T_avg, Y_g, q_g, d_tub, rel_rough, L, theta, M_g, R)


% Z factor block

P_pc = psc_press(Y_g);
T_pc = psc_tempt(Y_g);

P_pr = psr_press(P1, P_pc);
T_pr = psr_tempt(T_avg,T_pc);
Z_avg = Z_factor(T_pr, P_pr);

% formation volume block
B_g = gas_fom_vol(Z_avg, T_avg, P1);

% density block

rho_g = den_g(P1, M_g, Z_avg, R, T_avg); %lbm/ft3

% viscosity block
A = m_gA(M_g, T_avg);
B = m_gB(M_g, T_avg);
C = m_gC(B);
m_g = visz_gas(A, B, C, rho_g);

% Reynolds number block
if theta < 180
N_Re = Rey_num_ver(Y_g, q_g, d_tub, m_g);
else 
    N_Re = Rey_num_hor(Y_g, q_g, d_tub, m_g);
end

if N_Re <= 2100 %Laminar condition
    f_fact = 16/N_Re;
else %turbulent condition
    f_fact = friction_fact(N_Re, rel_rough);
end

if theta < 180
    s_P = s_P2(Y_g, L, theta, Z_avg, T_avg);
    P2 = P2_cal_ver(P1, s_P, f_fact, Z_avg, T_avg, q_g, theta, d_tub); %psi
else
    P2 = P2_cal_hor(P1, Y_g, f_fact, Z_avg, T_avg, q_g, d_tub);
end
end