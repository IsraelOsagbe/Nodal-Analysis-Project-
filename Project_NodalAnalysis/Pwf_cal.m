function Pwf = Pwf_cal(P_avg, q, Z, T, m_g, k, h, r_e, r_w)
Pwf = sqrt(P_avg^2 - (1424*q*Z*T*m_g*log(0.472*r_e/r_w)/k*h));
end