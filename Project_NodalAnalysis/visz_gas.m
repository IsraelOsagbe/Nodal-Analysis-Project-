%% 1 pound per cubic foot lb/ft3	= 0.016	grams per cubic centimeter g/cm3

function m_g = visz_gas(m_gA, m_gB, m_gC, rho_g)
m_g = m_gA.*(10^-4).*exp(m_gB.*((0.016.*rho_g).^m_gC)); %converting rho_g to g/cc
end