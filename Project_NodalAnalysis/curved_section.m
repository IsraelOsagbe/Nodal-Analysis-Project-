clear
close all
clc

main
Parameter

R = 3000;

MD_c = 1500*pi;
theta0 = 90;

sum_theta = 180;

theta1 = theta0 + ((360*L) / (2*pi*R));
theta = theta1- theta0;
Z1 = R*sind(theta);
j_theta = theta0:theta:sum_theta;


T_AHD = KOP:L:(KOP+MD_c); % from formular for length of a sector, theta = 90

P_cur = zeros(1,length(T_AHD));
T_cur = zeros(1,length(T_AHD));
x_cur = zeros(1,length(T_AHD));
y_cur = zeros(1,length(T_AHD));
Z_u = zeros(1,length(T_AHD));

P_avg_cur = zeros(1,length(T_AHD) - 1);
T_avg_cur = zeros(1,length(T_AHD) - 1);
T_pr_cur  = zeros(1,length(T_AHD) - 1);
P_pr_cur  = zeros(1,length(T_AHD) - 1);
Z_avg_cur = zeros(1,length(T_AHD) - 1);
B_g_cur = zeros(1, length(T_AHD) - 1);
rho_g_cur = zeros(1, length(T_AHD) - 1);
m_g_cur = zeros(1, length(T_AHD) - 1);



T_AHD(1) = KOP;
Z_u(1) = Z1;
j_theta(1) = theta;
T_cur(1) = T(21);
P_cur(1) = P(21);



r_MSE_cur = ones(1,length(T_AHD)-1);

for n = 2:length(T_AHD)
    for k = 2:length(j_theta) 
        Z_u(k) = (R*(sind(j_theta(k)- 90) - sind(j_theta(k-1) - 90))); %vertical height for deviating section
        y_cur(k) = KOP + R*(sind(j_theta(k) -90)); %vertical depth of deviated section
        y_cur(1) = KOP;
        x_cur(k) = R - (R*(cosd(j_theta(k) -90))); % horizontal distance from vertical point for deviated section
        x_cur(1) = 0;
    end
% end
% for k = 2:length(Z_u)
    T_cur(k) = T_grad*(Z_u(k-1) - Z_u(k)) + T_cur(k-1);
    T_avg_cur(k-1) = (T_cur(k-1) - T_cur(k))/(log(T_cur(k-1)/T_cur(k)));
    T_pr_cur(k-1) = psr_tempt(T_avg_cur(k-1),T_pc);

    while r_MSE_cur(k-1) > 0.005

        [P_cur(k), s_P_cur, f_fact_cur, N_Re_cur, m_g_cur(k), rho_g_cur(k), B_g_cur(k), Z_avg_cur(k)] = ultim_func(P_cur(k-1), T_avg_cur(k-1), Y_g, q_g, d_tub, rel_rough, L, j_theta(k-1), M_g, R);

        P_avg_cur(k-1) = (P_cur(k) + P_cur(k-1))/2;
        P_pr_cur(k-1)  = psr_press(P_avg_cur(k-1), P_pc);
        Z_avg_cur(k-1) = Z_factor(T_pr_cur(k-1), P_pr_cur(k-1));

        r_MSE_cur(k-1) = obj_func(Z_avg_cur(k-1), Z_avg(k));

    end
    B_g_cur(k-1) = gas_fom_vol(Z_avg_cur(k-1), T_avg_cur(k-1), P_avg_cur(k-1));
    rho_g_cur(k-1) = den_g(P_avg_cur(k-1), M_g, Z_avg_cur(k-1), R, T_avg_cur(k-1));
    m_g_cur(k-1) = visz_gas(m_gA(M_g, T_avg_cur(k-1)), m_gB(M_g, T_avg_cur(k-1)), m_gC(m_gB(M_g, T_avg_cur(k-1))), rho_g_cur(k-1));


end

plot (x_cur, y_cur, 'r-o')

