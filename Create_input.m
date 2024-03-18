clc
clear
close all
format long
addpath ( genpath ( 'Sample_specific_data/' ) ) 
addpath ( genpath ( 'Subroutines/Functions/' ) )
addpath ( genpath ( 'Subroutines/process_pdf' ) )


Name_of_input = 'First_example_324nm_0_17omega_14keV_to_18keV_Efilter';



% -------------------------------------------------------------------------
% INSTRUMENT
% -------------------------------------------------------------------------
instr = struct();

instr.goniometer_length = 240*1e7;    % [mm --> Å]
instr.wavelength        = 0.7107488;  % [Å]
instr.pol_angle         = 0;        % [deg]

instr.divergence        = 1/16;       % [deg]
instr.source_slit       = 0.1*1e7;   % [mm --> Å]
instr.det_slit          = 8.89*1e7;   % [mm --> Å]

instr.angle_in          = 0.15;      % Incidence angle [deg]

% -------------------------------------------------------------------------
% SAMPLE
% -------------------------------------------------------------------------
sample = struct();

sample.E   = ["V", "Zr"]; % Name of elements
sample.C   = [0.33, 0.67];      % Composition

sample.rho = 6.373;      % Density [g/cm^3]
sample.mu  = 1/(0.207204);  % Attenuation coefficient [1/microns]
sample.t   = 3240;       % Thickness [Å]
sample.len = 20*1e7;     % length sample [Å]


% -------------------------------------------------------------------------
% sub
% -------------------------------------------------------------------------

sample.sub.flag = ['y'];
sample.sub.E    = ["Al", "O"];
sample.sub.C    = [2, 3];

sample.air.E    = ["N", "O"];
sample.air.C    = [0.791, 0.209];

% -------------------------------------------------------------------------
% PDF process flags
% -------------------------------------------------------------------------
process = struct();

% XRD scan data path and loading procedure
process.load_scan_flag = 'process_pdf';  % ['XYE', 'process_pdf']
% On Windows, change the forward slash to a backslash
process.dir            = './Data_folder/First_example_324nm_0_17omega_14keV_to_18keV_Efilter/';

% background scan. Should be scaled by 0.5 I think. air_slit_1_16_new_energywindow_omega_0_13
% Background reduction
process.bkgr_flag = ["n", "process_pdf"];  % ['y', 'n']
process.bkgr_dir  = '';

% Apply geometric, polarization and absorption corrections
process.corr_flag = [0, 1, 1]; % [Geometric corr, Pol corr, Aborb corr]

% Dark counts
process.dark_I    = 0.0;

% Rebin flag
process.rebin_flag = ["n", 1];

% Lanczos filter
process.window_flag = 1;

% Dampening
process.damp_flag = 'n';
process.dampening = 0.0001;

% MC calculations type ['MC', 'spur_Q', 'spur_Q2]
% MC FOM type ["STD", "INT"]
process.MC_type_flag = ["MC", "INT", "incoh_VARY", "span_VARY"];

process.plot_flag    = ['y'];
% Define R
process.R = linspace(0.001, 14, 1000);

% -------------------------------------------------------------------------

[instr, sample] = set_up_input(instr, sample, process, Name_of_input);

function [instr, sample] = set_up_input(instr, sample, process, Name_of_input)


data = struct();


instr.source_beam_width = instr.source_slit + 0.1*instr.goniometer_length*2*tand(instr.divergence/2);

sample.C             = sample.C./sum(sample.C);
sample.mu            = sample.mu*1e-4;           %[1/microns -->1/Å]
[sample.Z, sample.M] = Read_elemental_data(sample.E);
sample.M_av          = sum(sample.M.*sample.C);
sample.FF_coeff      = Read_form_factor_coefficients(sample.Z, instr.wavelength);
sample.Compton_coeff = Read_Compton_coefficients(sample.Z);


[sample.air.Z, sample.air.M] = Read_elemental_data(sample.air.E);
sample.air.FF_coeff          = Read_form_factor_coefficients(sample.air.Z, instr.wavelength);

[sample.sub.Z, sample.sub.M] = Read_elemental_data(sample.sub.E);
sample.sub.Compton_coeff     = Read_Compton_coefficients(sample.sub.Z);
sample.sub.C                 = sample.sub.C./sum(sample.sub.C);
% -------------------------------------------------------------------------
% Load data
% -------------------------------------------------------------------------
if strcmp(process.load_scan_flag, 'process_pdf')
    disp('Loading data via the process_pdf method')
    [data.tth, data.I_scan, data.T,  data.I_error] = process_pdf(process.dir, instr.angle_in, 1);
elseif strcmp(process.load_scan_flag, 'process_pdf_cat')
    disp('Hello. This is your computer speaking. Please select another setting. Thank you.')
elseif strcmp(process.load_scan_flag, 'XYE')

     'ok'
     d = XRDMLread(strcat(process.dir, '.xrdml'));
     data.tth     = d.x';
     data.I_scan  = d.data';
     data.T       = d.time;
     data.I_error = sqrt(sqrt( d.data./d.time ))';
end

% Calculate Q range
data.Q = 4.*pi.*sind((data.tth)./2)./instr.wavelength;

% Calculate form factors
[data.f_sqrd, data.f_av_sqrd] = Form_factors(data.Q, sample.FF_coeff, sample.C);

[data.air.f_sqrd, data.air.f_av_sqrd] = Form_factors(data.Q, sample.air.FF_coeff, sample.air.C);
% Calcualte Compton scattering
data.I_compton = Compton_factors(data.Q, sample.Compton_coeff, sample.C, sample.Z, instr.wavelength);



data.sub.I_compton = Compton_factors(data.Q, sample.sub.Compton_coeff, sample.sub.C, sample.sub.Z, instr.wavelength);

% -------------------------------------------------------------------------
% Reduce background
% -------------------------------------------------------------------------
if strcmp(process.bkgr_flag(1), 'y')
    if strcmp(process.bkgr_flag(2), 'process_pdf')
        disp('Loading background via the process_pdf method')
        [data.tth_bkgr, data.I_bkgr, ~,  data.e_bkgr] = process_pdf(process.bkgr_dir, instr.angle_in, 1);
        data.I_red = data.I_scan - 0.5.*data.I_bkgr;
    else
        
    end
else
    data.I_red = data.I_scan;
end
data.I_red = data.I_red - process.dark_I.*data.T./data.T(1);

% -------------------------------------------------------------------------
% Apply corrections
% -------------------------------------------------------------------------
if process.corr_flag(1) == 1
    disp('fix')
    %data.corr.geo = Geometric_factor_correction(instr.angle_in, data.tth - instr.angle_in, instr.source_beam_width, instr.det_slit);
%     data.corr.geo =  sin((data.tth + 90)./90).^0.5;
    %data.corr.geo =  sin((data.Q + 4*pi)./(4*pi));
    %data.corr.geo =  sind(data.tth./2 + data.tth(end)/2);
    %data.corr.geo =  sind(data.tth + 60);
    

%     frac = ones(length(data.tth), 1);
%     rot = 0.1;
%     L              = sample.len/2;
%     det_foot       = instr.det_slit./sind(data.tth./2);
%     rot_sample_len = (L*(sqrt(2)*cosd(45 - rot) - 1) + L);
%     foot_ind = det_foot <  rot_sample_len;
% 
%     D = sqrt(2).*cosd(45 - rot) - det_foot(foot_ind)./L;
%   
%     
%     tan_fact = 1;% tand(rot) + tand(90 - rot);
%     frac(foot_ind) = 1 -D.^2.*tan_fact./(1 - D.^2.*tan_fact);


    %sample_footprint = (instr.source_slit)/sind(instr.angle_in)/2;
    %geo_frac  = instr.det_slit/(sample_footprint);
%     geo_frac  = instr.det_slit/(sample.len/2);
%     angle_geo = asind(geo_frac);
%     angle_ind = data.tth - instr.angle_in >= angle_geo;
%     data.corr.geo = ones(length(data.tth), 1);
%     data.corr.geo(angle_ind) = geo_frac./sind(data.tth(angle_ind) - instr.angle_in);

    sample_footprint = (instr.source_slit)/sind(instr.angle_in)/2;
    geo_frac  = instr.det_slit/(sample_footprint/2);

    bw = 0.02637*1e7;
    bw = 0.1*1e7;

    F = bw/sind(instr.angle_in);
    F(F > sample.len) = sample.len;

    Fd = instr.det_slit./sind(data.tth - instr.angle_in);
    Fd(Fd>sample.len) = sample.len;

%     data.corr.geo = Fd./F;

    geo_frac  = instr.det_slit./F
    angle_geo = asind(geo_frac)
    angle_ind = data.tth - instr.angle_in >= angle_geo;
    data.corr.geo = ones(length(data.tth), 1);
    data.corr.geo(angle_ind) = geo_frac./sind(data.tth(angle_ind) - instr.angle_in);

%     figure
%     hold on
%     plot(data.tth, data.corr.geo, 'k')
%     plot(data.tth, frac, 'r')
% 
%     asdsadsadd
else
    data.corr.geo = ones(length(data.tth), 1);
end

if process.corr_flag(2) == 1
    data.corr.pol = Polarisation_factor_correction(data.tth, instr.pol_angle);
else
    data.corr.pol = ones(length(data.tth), 1);
end

if process.corr_flag(3) == 1
    data.corr.abs = Num_of_scatters(instr.angle_in, data.tth, sample.t, sample.mu);
else
    data.corr.abs = ones(length(data.tth), 1);
end

data.corr.tot = data.corr.geo.*data.corr.pol.*data.corr.abs;

data.I_corr       = data.I_red./data.corr.tot;
data.I_corr_error = data.I_error./data.corr.tot;

% -------------------------------------------------------------------------
if strcmp(process.bkgr_flag(1), 'y')
    plot_I_bkgr(data)
end
save(strcat('Sample_specific_data/Input_data/', Name_of_input),  'instr', 'sample', 'process', 'data')
disp('Done')
% 
% data.I_corr = data.I_corr.*data.corr.geo;
% 
%     geo_frac  = instr.det_slit/(1.*sample.len);
%     angle_geo = asind(geo_frac);
%     angle_ind = data.tth - instr.angle_in >= angle_geo;
%     data.corr.geo = ones(length(data.tth), 1);
%     data.corr.geo(angle_ind) = geo_frac./sind(data.tth(angle_ind) - instr.angle_in);
% 
% data.I_corr = data.I_corr./data.corr.geo;
% 
% figure
% hold on
% plot(data.tth, data.I_corr)
% 


% asdsad
end


function plot_I_bkgr(data)
figure
hold on
plot(data.tth, data.I_scan, 'k')
plot(data.tth_bkgr, data.I_bkgr, 'r')
plot(data.tth, data.I_red, 'b')
end