
function [rho, g] = Fit_sine(R, G, R_guess)



[G_slope, min_i] = Find_peak_minimum(R, G, R_guess);

R_osc = transpose(R(1:min_i));
G_osc = G(1:min_i) + R_osc;

%FT = fittype('a.*sin(b.*x) + c.*cos(b.*x)');
%FT = fittype('a.*sin(b.*x) + c.*x');
%FT = fittype('a.*sin(b.*x).*exp(-c.*x)');
FT = fittype('a.*sin(b.*x)');


[fit_data, ~, ~] = fit(R_osc, G_osc, FT, 'StartPoint', [1,12]);
K = coeffvalues(fit_data)
%sine_data = K(1).*sin(K(2).*R_osc) + K(3).*R_osc;
%sine_data = K(1).*sin(K(2).*R_osc).*exp(-K(3).*R_osc);
sine_data = K(1).*sin(K(2).*R_osc);
figure
hold on
%plot(R_osc, G_osc, 'k')
%plot(R_osc, sine_data, 'r')

sine_data = K(1).*sin(K(2).*R).*exp(-0.25.*R_guess.*R);
 plot(R, G, 'k')
 plot(R, sine_data, 'r')
 plot(R, G - sine_data', 'b--')

end