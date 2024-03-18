
function [N_coh_ref, N_inc_ref, std_ref] = Monte_Carlo_nomalization_sum_rule(Q, I, f_sqrd, f_av_sqrd, I_compton, num_iterations, span_interval, guess_of_first_peak,...
                                                                    N_coh_ref, N_inc_ref, window_flag, damp_flag, damping, M_av, rho_ref)
fprintf('\n')
fprintf([' Q max: ', num2str(Q(end)), ' \n'])

N_coh   = N_coh_ref;
N_inc   = N_inc_ref;
std_ref = 5;
N_A     = 6.022*1e23;

R = linspace(0.001, guess_of_first_peak + 1, 50);

figure
hold on
for i = 1:num_iterations
    
    [~, ~, F, ~, G]  = Intensity_to_F_to_PDF_calc(Q, R, I, f_sqrd, f_av_sqrd, I_compton, N_coh, N_inc, window_flag, damp_flag, damping);
    
    
    
    rho0 = sum_rule(Q, F);
    rho  = rho0*M_av/N_A*1e24;
    std_current      = abs(rho - rho_ref);
    
    
    % If better std, store the new parameters and decrease the search
    % interval 
    if std_current < std_ref
        std_ref        = std_current;
        N_coh_ref      = N_coh;
        N_inc_ref      = N_inc;
        %span_interval  = span_interval/2;
        plot(R, G)
    end
    
    N_coh = abs( span_interval*(rand - 0.5) + N_coh_ref );
    N_inc = abs( span_interval*(rand - 0.5) + N_inc_ref );
    
    if mod(i, round(num_iterations/10)) == 0
        fprintf([num2str(100*i/num_iterations), ' percent done \n'])
    end
    
end
xlabel('R')
ylabel('G(R)')
hold off
pause(0.5)
end



