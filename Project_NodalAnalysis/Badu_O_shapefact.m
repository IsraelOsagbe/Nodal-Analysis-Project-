function ln_CH = Badu_O_shapefact(a,h,I_ani,y0_a,z0_h)

ln_CH = (6.28*(a/(I_ani*h))*((1/3)-(y0_a)+((y0_a)^2))) - log(sin(pi*z0_h)) -0.5*(a/(I_ani*h)) - 1.088;

end