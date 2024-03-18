clc
clear
close all
format long
tic

addpath ( genpath ( 'Theory_SIO2/' ) )
addpath ( genpath ( 'Subroutines/' ) )
addpath ( genpath ( 'Sample_specific_data/Input_data' ) )
MC = struct();

% name_input should be the name used in Create Input
name_input = 'First_example_324nm_0_17omega_14keV_to_18keV_Efilter';
%name_input = 'Second_example_324nm_0_17omega_17_15keV_to_18keV_Efilter';

load(name_input)

MC_input      = load('MC_input.txt');
MC.num_steps  = MC_input(1);
MC.span(1)    = MC_input(2);
MC.span(2)    = MC_input(3);
MC.span(3)    = MC_input(4);
MC.peak_guess = MC_input(5);
MC.Q_max      = [MC_input(6):MC_input(8):MC_input(7)];
N_coh_best    = MC_input(9);
N_inc_best(1) = MC_input(10);
N_inc_best(2) = MC_input(11);

MC.span(4)    = MC_input(12);
MC.span(5)    = MC_input(13);
N_inc_best(3) = MC_input(14);
N_inc_best(4) = MC_input(15);

if strcmp(process.MC_type_flag(3), "incoh_FIXED")
    N_inc_best(1) = 1;
    N_inc_best(2) = 1;
    N_inc_best(3) = 1;
end


process.window_flag = 0;



Best_value = struct('std', 1000, 'Q_best', 0, 'N_coh_best', N_coh_best, 'N_inc_best', N_inc_best, 'vec',[0, 0, 0, N_coh_best, N_inc_best]);


for index = 1:length(MC.Q_max)
    Q_max            = MC.Q_max(index);
    [~, max_guess_i] = min(abs(data.Q - Q_max));
    I_corr           = data.I_corr(1:max_guess_i);
    Q                = data.Q(1:max_guess_i);
    f_sqrd           = data.f_sqrd(1:max_guess_i);
    f_av_sqrd        = data.f_av_sqrd(1:max_guess_i);
    I_compton        = data.I_compton(1:max_guess_i);
    if strcmp(sample.sub.flag(1), 'y')
        I_compton_sub    = data.sub.I_compton(1:max_guess_i);
        f_sqrd_air       = data.air.f_sqrd(1:max_guess_i);
        f_av_sqrd_air       = data.air.f_av_sqrd(1:max_guess_i);

        f_sqrd_air  = f_sqrd_air./f_av_sqrd_air;
    end

    % Calculate the N_coh and N_inc via randomly selecting values within a span
    % around the guess value. A good start is to have a very large span and
    % after several iterations use a small interval to hone in to the better
    % precicise the values
    % Parameters that goes into function:
    % Q, I_corrected, f_sqrd, f_av_sqrd, I_compton, num_iterations, Span_of_search_interval, guess_of_first_peak, N_coh_best, N_inc_best)

    if strcmp(sample.sub.flag(1), 'y')
%         [N_coh, N_inc, std_ref] = Monte_Carlo_nomalization_sub(Q, I_corr, f_sqrd, f_av_sqrd, I_compton, I_compton_sub, MC.num_steps, MC.span, MC.peak_guess, ...
%             Best_value.N_coh_best, Best_value.N_inc_best, process);

        [N_coh, N_inc, std_ref] = Monte_Carlo_nomalization_sub_EfilterVary(Q, I_corr, f_sqrd, f_av_sqrd, I_compton, I_compton_sub, MC.num_steps, MC.span, MC.peak_guess, ...
            Best_value.N_coh_best, Best_value.N_inc_best, process, instr.wavelength);

        %     [N_coh, N_inc, std_ref] = Monte_Carlo_nomalization_sub_air(Q, I_corr, f_sqrd, f_av_sqrd, f_sqrd_air, I_compton, I_compton_sub, MC.num_steps, MC.span, MC.peak_guess, ...
        %         Best_value.N_coh_best, Best_value.N_inc_best, process);
    else
        [N_coh, N_inc, std_ref] = Monte_Carlo_nomalization(Q, I_corr, f_sqrd, f_av_sqrd, I_compton, MC.num_steps, MC.span, MC.peak_guess, ...
            Best_value.N_coh_best, Best_value.N_inc_best, process);
    end

    if std_ref < Best_value.std
        Best_value.std        = std_ref;
        Best_value.Q_best     = Q_max;
        Best_value.N_coh_best = N_coh;
        Best_value.N_inc_best = N_inc;
    end
    Best_value.vec = [Best_value.vec; [index, std_ref, Q_max, N_coh, N_inc]];

end

data.N_coh = Best_value.N_coh_best;
data.N_inc = Best_value.N_inc_best;
[~, max_guess_i] = min(abs(data.Q - Best_value.Q_best));
data.I_corr_orig = data.I_corr;
data.I_corr      = data.I_corr(1:max_guess_i);
data.Q           = data.Q(1:max_guess_i);
data.f_sqrd      = data.f_sqrd(1:max_guess_i);
data.f_av_sqrd   = data.f_av_sqrd(1:max_guess_i);
data.I_compton   = data.I_compton(1:max_guess_i);

if strcmp(sample.sub.flag(1), 'y')
    data.sub.I_compton   = data.sub.I_compton(1:max_guess_i);
    data.air.f_sqrd      = data.air.f_sqrd(1:max_guess_i);
    data.air.f_av_sqrd   = data.air.f_av_sqrd(1:max_guess_i);

    data.air.f_sqrd      = data.air.f_sqrd./data.air.f_av_sqrd;
end
data.R           = process.R;

% From the calculated N_coh and N_inc, calculate a final structure factor
% and reduced PDF.
% Parameters that goes into function:
% Q, R, I_corrected, f_sqrd, f_av_sqrd, I_compton, N_coh, N_inc, Damping)

if strcmp(sample.sub.flag(1), 'y')

    process.window_flag = 0;
%     [data.Q_inter, data.I_norm, data.F_orig, data.F, data.G] = ...
%         Intensity_to_F_to_PDF_calc_sub(data.Q, data.R, data.I_corr, data.f_sqrd, data.f_av_sqrd, data.I_compton, data.sub.I_compton, data.N_coh, data.N_inc, process);
    [data.Q_inter, data.I_norm, data.F_orig, data.F, data.G] = ...
    Intensity_to_F_to_PDF_calc_sub_EfilterVary(data.Q, data.R, data.I_corr, data.f_sqrd, data.f_av_sqrd, data.I_compton, data.sub.I_compton, data.N_coh, data.N_inc, process, instr.wavelength);

    process.window_flag = 1;
%     [~, ~, ~, ~, data.G_L] = ...
%         Intensity_to_F_to_PDF_calc_sub(data.Q, data.R, data.I_corr, data.f_sqrd, data.f_av_sqrd, data.I_compton, data.sub.I_compton, data.N_coh, data.N_inc, process);
    [~, ~, ~, ~, data.G_L] = ...
    Intensity_to_F_to_PDF_calc_sub_EfilterVary(data.Q, data.R, data.I_corr, data.f_sqrd, data.f_av_sqrd, data.I_compton, data.sub.I_compton, data.N_coh, data.N_inc, process, instr.wavelength);
    % [data.Q_inter, data.I_norm, data.F_orig, data.F, data.G] = ...
    %     Intensity_to_F_to_PDF_calc_sub_air(data.Q, data.R, data.I_corr, data.f_sqrd, data.f_av_sqrd, data.air.f_sqrd, data.I_compton, data.sub.I_compton, data.N_coh, data.N_inc, process);
else
    [data.Q_inter, data.I_norm, data.F_orig, data.F, data.G] = ...
        Intensity_to_F_to_PDF_calc(data.Q, data.R, data.I_corr, data.f_sqrd, data.f_av_sqrd, data.I_compton, data.N_coh, data.N_inc, process);
end

if data.N_inc(3) > 0
%         angle_EFilter = asind(data.N_inc(3));
%         angle_ind = data.tth >= angle_EFilter;
%         data.corr.geo = ones(length(data.tth), 1);
%         data.corr.geo(angle_ind) = data.N_inc(3)./sind(data.tth(angle_ind));


    Q_shift = data.N_inc(4).*max(data.Q);
    data.corr.geo = 1./(exp( data.N_inc(3).*(data.Q - Q_shift)) + 1) ;
end
% Error
data.I_norm_error = data.N_coh.*data.I_corr_error(1:max_guess_i);
data.F_error      = data.Q.*data.I_norm_error./data.f_av_sqrd;

% Convert G --> g and calculate the density via the slope
% Parameters that goes into function:
% R, G, guess_peak_position
[data.rho, data.g] = Calc_g_and_density(data.R, data.G, MC.peak_guess);
[data.rho, data.g_L] = Calc_g_and_density(data.R, data.G_L, MC.peak_guess);

% Calculate the coordination number from a selected peak
% Parameters that goes into function:
% R, G, guess_peak_position
N_A = 6.022*1e23;
density = sample.rho/(sample.M_av/N_A*1e24);
%Coordination_num = Calc_coordination_number(data.R, data.g, data.rho, MC.peak_guess)



%rho0 = sum_rule(data.Q_inter, data.F);
rho0 = sum_rule(data.Q, data.F_orig);

disp(['Sum rule gives density: ',num2str(rho0*sample.M_av/N_A*1e24),' g/cm^3'])

% density_layers_dummy*N_A/mole_mass_dummy
data.rho_slope = data.rho*sample.M_av/N_A*1e24; % g/cm^3
disp(['Slope of g(r) gives density: ',num2str(data.rho_slope),' g/cm^3'])


data.N_coh_standard = Normalization_warren(data.Q, data.I_corr, data.f_sqrd, data.f_av_sqrd, data.I_compton, sample.rho, data.N_inc, process.MC_type_flag);
[data.beta, data.beta_slope, data.beta_warren, data.beta_vec] = Calc_beta(instr, data, sample, Best_value);

disp([' '])
disp(['Beta: ',num2str(data.beta, 6), ''])
disp(['Beta Norm: ',num2str(data.beta_warren, 6), ''])
disp(['Beta from slope density: ',num2str(data.beta_slope, 6), ''])
disp([' '])
% disp(['Beta from Theory/Calc_beta_siO2: ',num2str(1.920828733250473, 6),'e+15 +-', num2str(0.04755899213181691, 6), '+e15 '])
% disp(['Window Beta from Theory/Calc_beta_siO2: ',num2str(1.895780708534589, 6),'e+15 +-', num2str(0.03970920890367402, 6), '+e15'])
%disp(['Beta Warren: ',num2str(1.66532854754341, 6), 'e+15'])
disp(['Beta Warren: ',num2str(3.008517258131464, 6), 'e+12'])

run Subroutines/Plot_graphs

fid   = fopen('MC_input.txt');
file  = textscan(fid,'%s','delimiter','\n');
MC.file = file{1};
fclose(fid);
MC.file{14} = num2str(Best_value.Q_best);
MC.file{15} = num2str(Best_value.Q_best);
MC.file{19} = num2str(data.N_coh);
MC.file{22} = num2str(data.N_inc(1));
MC.file{23} = num2str(data.N_inc(2));


MC.file{30} = num2str(data.N_inc(3));
MC.file{31} = num2str(data.N_inc(4));

writecell(MC.file, 'MC_input.txt', "QuoteStrings", 0)
toc

