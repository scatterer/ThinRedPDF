
function [Q, Q_inter, I_norm, F_orig, F, G] = Intensity_to_F_to_PDF_calc(Q, R, I, f_sqrd, f_av_sqrd, I_compton, N_coh, N_inc, window_flag, damp_flag, alpha)


I_norm         = N_coh.*I - N_inc.*I_compton;
F_orig         = (I_norm - f_sqrd)./f_av_sqrd;
[Q_inter, F]   = F_interpolate(Q, F_orig);


Q_new          = linspace(Q_inter(1), Q_inter(end), length(Q_inter)/2);
[Q_new, F_new] = rebin(Q_new, Q_inter, F);
Q_inter        = Q_new;

if window_flag == 1
    F = F_new.*sin(pi.*Q_inter./Q_inter(end))./(pi.*Q_inter./Q_inter(end));
else
    F = F_new;
end

if damp_flag == 1
    G = PDF_calculation(Q_inter, R, F.*exp(-alpha.^2.*Q_inter.^2));
else
    G = PDF_calculation(Q_inter, R, F);
end
end


