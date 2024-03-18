
function geo_corr_factor = Geometric_factor_correction(angle_in, angle_out, source_beam_width, det_slit_width)
geo_corr_factor = zeros(length(angle_out), 1);


g_factor = (det_slit_width/source_beam_width).*(sind(angle_in)./sind(angle_out));

for i = 1:length(angle_out)
    if g_factor(i) < 1
        geo_corr_factor(i) = g_factor(i);
    else
        geo_corr_factor(i) = 1;
    end
end

end