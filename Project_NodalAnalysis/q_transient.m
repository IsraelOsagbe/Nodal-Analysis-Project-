function q = q_transient(k, h, P_avg, Pwf, poro, T, t, c_g, r_w, m_g)
q = (k*h*(P_avg^2 - Pwf^2))/(1638*T* (log10(t) + log(k/(poro*m_g*c_g*r_w^2))) - 3.23);
end