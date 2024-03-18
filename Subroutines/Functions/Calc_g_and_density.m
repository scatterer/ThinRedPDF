
function [rho, g] = Calc_g_and_density(R, G, R_guess)



[G_slope, min_i] = Find_peak_minimum(R, G, R_guess);

rho = abs( -G_slope(end)./(4.*pi.*R(min_i)) ) ;


g = G./(4*pi.*rho.*transpose(R)) + 1;

end