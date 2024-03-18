
function num_scatter = Num_of_scatters(angle_in, angle_out, thickness, mu)

coeff               = (1 + sind(angle_in)./sind(angle_out));
exponent_angle_part = 1./sind(angle_in) + 1./sind(angle_out);
beer_lambert_part   = 1 - exp(-mu.*thickness.*exponent_angle_part);
num_scatter         = (1./coeff).*beer_lambert_part;

end