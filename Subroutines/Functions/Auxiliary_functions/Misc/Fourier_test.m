% [Q, F, S] = Calc_Structure_Factor_Fourier(data.R', data.G, data.Q_inter, 0, MC.peak_guess);
% [Q, F_red, S] = Calc_Structure_Factor_Fourier(data.R', data.G, data.Q_inter, 1, MC.peak_guess);
% 
% F_sub =  F - F_red;
% F_sub(1) = F_sub(2);
% G = PDF_calculation(Q, data.R, F_sub);
% [Q, F, S] = Calc_Structure_Factor_Fourier(data.R', G, Q, 0, MC.peak_guess);
% [Q, F_red, S] = Calc_Structure_Factor_Fourier(data.R', G, Q, 1, MC.peak_guess);
% F_sub =  F - F_red;
% F_sub(1) = F_sub(2);
% G = PDF_calculation(Q, data.R, F_sub);
% 
% for x = 1:1
% [Q, F, S] = Calc_Structure_Factor_Fourier(data.R', G, Q, 0, MC.peak_guess);
% [Q, F_red, S] = Calc_Structure_Factor_Fourier(data.R', G, Q, 1, MC.peak_guess);
% F_sub =  F - F_red;
% F_sub(1) = F_sub(2);
% G = PDF_calculation(Q, data.R, F_sub);
% end
% 
% data.G = G;
% data.F = F_sub;
% 
% figure
% hold on
% plot(data.Q_inter, data.F.*data.Q_inter,'-')
% plot(Q, F.*Q,'-')
% plot(Q, F_red.*Q,'-')
% 
% 
% figure
% hold on
% plot(data.R, data.G, 'o')
% plot(data.R, G,'o')
% xlim([0,6])