clc
clear
close all
format long

addpath ( genpath ( 'Subroutines/Functions/' ) )
addpath ( genpath ( 'Subroutines/' ) )
run Sample_specific_data/D8_data
%run Sample_specific_data/SiO2 

run Sample_specific_data/H7_V50_Zr50

% Path and name of scan file and background

%file_name_scan = 'Diffraction_data/SiO2_005mmts_nonifilter/gixrd_exported.xye';
%file_name_bkgr = 'Diffraction_data/SiO2_005mmts_nonifilter/air_gi_exported.xye';

file_name_scan = 'Diffraction_data/PDF_210205/gixrd_full_detector_0D_0_828deg_vct_450s_exported.xye';
file_name_bkgr = 'Diffraction_data/PDF_210205/bkg_0_828_exported.xye';

% load data
run Subroutines/Read_scan_data

% Background reduction and angle dependent correction factors
%I_red               = bkgr_reduction(theta_scan, I_scan, theta_bkgr, I_bkgr);
angle_out           = theta_scan - angle_in;


I_reduced = substrate_reduction_func(I_scan, I_bkgr, angle_in, angle_out, mu, thickness);

figure
plot(theta_scan, I_scan)
hold on
plot(theta_scan, I_reduced)

function y_reduced = substrate_reduction_func(y_scan, y_substrate, angle_in, angle_out, mu, thickness)



exponent_angle_part = 1./sind(angle_in) + 1./sind(angle_out);
beer_lambert_part   = exp(-mu.*thickness.*exponent_angle_part);
y_reduced           = y_scan - y_substrate.*beer_lambert_part;

end