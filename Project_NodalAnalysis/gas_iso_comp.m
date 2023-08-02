function c_g = gas_iso_comp(Z, P, Z_grad)
c_g = 1/P - (1/Z)*(Z_grad);
end