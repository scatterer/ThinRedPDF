
function pol_corr_factor = Polarisation_factor_correction(tth, pol_angle)

% monocromator source side 
pol_corr_factor = (1 + cosd(tth).^2.*cosd(2*pol_angle).^2 )./(1 + cosd(2*pol_angle).^2);

%pol_corr_factor = 0.5.*(1 + cosd(tth).^2.*cosd(2*pol_angle).^2);

end