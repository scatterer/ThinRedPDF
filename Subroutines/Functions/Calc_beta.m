function [beta, beta_slope, beta_standard, beta_vec] = Calc_beta(instr, data, sample, Best_value)

N_A   = 6.022*1e23;



t     = 1/sample.mu;
rho_0 = sample.rho/(sample.M_av/N_A*1e24);
rho_slope = data.rho_slope/(sample.M_av/N_A*1e24);

% from \Theory_SIO2\Calc_beta_SiO2
f = 0.1;

instr.source_beam_width = instr.source_slit + f*instr.goniometer_length*2*tand(instr.divergence/2);
%instr.source_beam_width = instr.goniometer_length*2*tand(instr.divergence/2);
sigma = instr.source_beam_width*sample.len;
gamma   = 0.28*pi/180; %instr.det_slit;


if length(Best_value.vec(:,1)) > 2
    N_coh    = Best_value.vec(2:end,4);
    beta_vec = N_coh*gamma*sigma*rho_0*t;
else
    beta_vec = 0;
end
N_coh = data.N_coh;
beta = N_coh*gamma*sigma*rho_0*t;
beta_slope = N_coh*gamma*sigma*rho_slope*t;
beta_standard = data.N_coh_standard*gamma*sigma*rho_0*t;
%beta = N_coh*sigma*rho_0*t

% warren = load('Theory_SIO2/Warren_Structure_factor_normalized.xy');
%  Q_warren = warren(1:end,1);
%  F_warren = warren(1:end,2)./Q_warren./2;
%  [f_sqrd, f_av_sqrd] = Form_factors(Q_warren, sample.FF_coeff, sample.C);
% 
% % Calcualte Compton scattering
% I_compton = Compton_factors(Q_warren, sample.Compton_coeff, sample.C, sample.Z);
% 
% % F_orig         = (I_norm - f_sqrd)./f_av_sqrd;
% I_coh_Warren = F_warren.*f_av_sqrd + f_sqrd  + 1.*I_compton;
% 
% [~,min_i]   =  min(abs(data.Q - 4));
% [~,min_i_W] =  min(abs(Q_warren- 4));
% 
% [~,max_i]   =  max(abs(data.Q));
% [~,max_i_W] =  min(abs(Q_warren- data.Q(max_i)));
% 
% N_Warren    = mean(I_coh_Warren(min_i_W:max_i_W))./mean(data.I_corr(min_i:max_i))
% beta_Warren = N_Warren*gamma*sigma*rho_0*t
% beta_gamma_Warren = N_Warren*sigma*rho_0*t
% 
% 
% [Q_warren_rebin, I_coh_Warren_rebin]= rebin(data.Q, Q_warren(1:max_i_W+1), I_coh_Warren(1:max_i_W+1));
% figure
% plot(Q_warren, I_coh_Warren,'k')
% hold on
% plot(data.Q, N_Warren.*data.I_corr,'r')
% 
% 
% ratio = N_Warren.*data.I_corr(1:end-1)./I_coh_Warren_rebin;
% 
% mean(ratio(min_i:end))
% figure
% plot(Q_warren_rebin(min_i:end), ratio(min_i:end),'k')




end
