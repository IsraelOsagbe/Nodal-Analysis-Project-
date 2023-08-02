function J = Badu_Odeh (k_y, k_z, b, m_o, b_o, cross_area, r_w, s, s_R, ln_CH)

J = (sqrt(k_y*k_z)*b)/(141.2*m_o*b_o*((log((cross_area^0.5)/r_w)) + log(ln_CH) - 0.75 + s_R + s));

end