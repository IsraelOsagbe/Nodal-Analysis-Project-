function q = q_cal (P1, P2, h, Z_avg, T_avg, m_g_avg, r_e, r_w, k_h)
q = ((k_h*h*(P2^2 - P1^2))/(1424*Z_avg*T_avg*m_g_avg*(log((0.472*r_e)/r_w))));
end