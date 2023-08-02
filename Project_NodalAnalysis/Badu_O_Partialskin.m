function s_R = Badu_O_Partialskin(v1, v2, v3, P_xyz, P_y, P_xy, P_xy_prime)

if (v1>v2)&&(v2>v3)
    s_R = P_xyz + P_xy_prime;
else
    
    s_R = P_xyz + P_y + P_xy;
end
end