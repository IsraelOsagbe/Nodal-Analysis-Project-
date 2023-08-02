function Ff = friction_fact(N_Re, rel_rough)
Ff = (1./(-4*log10((rel_rough/3.7065) - ((5.0452./N_Re).*log10(rel_rough^1.1098/2.8257 + ((7.149./N_Re).^0.8981)))))).^2;
end