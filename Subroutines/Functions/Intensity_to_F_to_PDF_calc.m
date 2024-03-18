
function [Q_inter, I_norm, F_orig, F, G] = Intensity_to_F_to_PDF_calc(Q, R, I, f_sqrd, f_av_sqrd, I_compton, N_coh, N_inc, process)

if strcmp(process.MC_type_flag(1), "spur_Q2")
    I_norm         = N_coh.*I - 1.*I_compton - N_inc.*Q.^2;
    %I_norm         = N_coh.*(I - 1.*I_compton) - N_inc.*Q.^2;
elseif strcmp(process.MC_type_flag(1), "spur_Q")
    I_norm         = N_coh.*I - 1*I_compton - N_inc.*Q;
    %I_norm         = N_coh.*(I - 1.*I_compton) - N_inc.*Q;
    
elseif strcmp(process.MC_type_flag(1), "spur_sine")
    %geo = sin((Q + N_inc.*(4*pi))./(N_inc.*(4*pi)));
    %geo = sin((Q + (4*pi))./((4*pi))).^N_inc;
    th = asind(Q./(4*pi)*0.7107488);
    geo = sind(th + N_inc);
    I_norm         = N_coh.*I./geo - 3.7325*I_compton;
    
elseif strcmp(process.MC_type_flag(1), "spur_pol")
    tth = 2.*asind(Q./(4*pi)*0.7107488);
    pol = Polarisation_factor_correction(tth, N_inc);
    I_norm         = N_coh.*I./pol - 1*I_compton;
    
else
    I_norm         = N_coh.*I - N_inc(1).*I_compton;
    %I_norm         = N_coh.*(I - N_inc.*I_compton);
end
F_orig         = (I_norm - f_sqrd)./f_av_sqrd;
[Q_inter, F]   = F_interpolate(Q, F_orig);

if strcmp(process.rebin_flag, 'y')
    Q_new          = linspace(Q_inter(1), Q_inter(end), length(Q_inter)/process.rebin_flag(2));
    [Q_new, F_new] = rebin(Q_new, Q_inter, F);
    Q_inter        = Q_new;
    F              = F_new;
end



if process.window_flag == 1
    F = F.*sin(pi.*Q_inter./Q_inter(end))./(pi.*Q_inter./Q_inter(end));
end

%warren = load('Theory_SIO2/Warren_Structure_factor_normalized.xy');
 %Q_inter = warren(1:2700,1);
 %Q_inter = warren(1:1920,1);
 %Q_inter = warren(1:2450,1);
 %F = warren(1:2700,2)./Q_inter./2;
 %F = warren(1:2450,2)./Q_inter./2;
%   Q_inter = warren(1:end,1);
%  F = warren(1:end,2)./Q_inter./2.*sin(pi.*Q_inter./Q_inter(end))./(pi.*Q_inter./Q_inter(end));

if strcmp(process.damp_flag, 'y')
    G = PDF_calculation(Q_inter, R, F.*exp(-process.dampening.^2.*Q_inter.^2));
else
    G = PDF_calculation(Q_inter, R, F);
end
end


