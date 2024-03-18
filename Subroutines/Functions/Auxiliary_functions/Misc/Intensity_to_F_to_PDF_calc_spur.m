
function [Q_inter, I_norm, F_orig, F, G] = Intensity_to_F_to_PDF_calc_spur(Q, R, I, f_sqrd, f_av_sqrd, I_compton, N_coh, N_inc, N_spur, window_flag, damp_flag, alpha)


I_norm         = N_coh.*I - N_inc.*I_compton - N_spur.*Q.^2;
F_orig         = (I_norm - f_sqrd)./f_av_sqrd;
[Q_inter, F]   = F_interpolate(Q, F_orig);




% smear_width = 0.025;
% F_new = transpose(GaussConv(Q_inter, F, smear_width*2.34));
% F_new(1:10) =  F(1:10);
% Q_new       = Q_inter(1):0.01:Q_inter(end);
% F_new_inter = interp1(Q_inter,F_new,Q_new);
% Q_inter        = Q_new;
% F_new = F_new_inter;

Q_new          = linspace(Q_inter(1), Q_inter(end), length(Q_inter)/2);
[Q_new, F_new] = rebin(Q_new, Q_inter, F);
Q_inter        = Q_new;

if window_flag == 1
    F = F_new.*sin(pi.*Q_inter./Q_inter(end))./(pi.*Q_inter./Q_inter(end));
else
    F = F_new;
end

if strcmp(damp_flag, 'y')
    G = PDF_calculation(Q_inter, R, F.*exp(-alpha.^2.*Q_inter.^2));
else
    G = PDF_calculation(Q_inter, R, F);
end
end


