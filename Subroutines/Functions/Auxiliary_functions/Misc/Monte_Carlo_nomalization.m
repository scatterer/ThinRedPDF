
function [N_coh_ref, N_inc_ref, std_ref] = Monte_Carlo_nomalization(Q, I, f_sqrd, f_av_sqrd, I_compton, num_iterations, span_interval, guess_of_first_peak,...
                                                                    N_coh_ref, N_inc_ref, window_flag, damp_flag, damping)
fprintf('\n')
fprintf([' Q max: ', num2str(Q(end)), ' \n'])

N_coh     = N_coh_ref;
N_inc     = N_inc_ref;
std_ref   = 5;

R = linspace(0.001, guess_of_first_peak + 1, 100);

figure
hold on
for i = 1:num_iterations
    
    [~, ~, ~, ~, G]  = Intensity_to_F_to_PDF_calc(Q, R, I, f_sqrd, f_av_sqrd, I_compton, N_coh, N_inc, window_flag, damp_flag, damping);
    
    [G_slope, min_i] = Find_peak_minimum(R, G, guess_of_first_peak);
    
    std_current      = std(G_slope - G(1:min_i));
    
    
    % If better std, store the new parameters and decrease the search
    % interval 
    if std_current < std_ref
        std_ref        = std_current;
        N_coh_ref      = N_coh;
        N_inc_ref      = N_inc;
        %span_interval  = span_interval/2;
        plot(R, G)
        plot(R(1:min_i), G_slope, '--k')
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
pause(1)
end



