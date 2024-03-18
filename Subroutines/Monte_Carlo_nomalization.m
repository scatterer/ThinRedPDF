
function [N_coh_ref, N_inc_ref, std_ref] = Monte_Carlo_nomalization(Q, I, f_sqrd, f_av_sqrd, I_compton, num_iterations, span_interval_ref, guess_of_first_peak,...
    N_coh_ref, N_inc_ref, process)
fprintf('\n')
fprintf([' Q max: ', num2str(Q(end)), ' \n'])

R             = linspace(0.001, guess_of_first_peak + 1, 100);
std_ref       = 5;
N_coh         = N_coh_ref;
N_inc         = N_inc_ref;
span_interval = span_interval_ref;

if strcmp(process.plot_flag(1), "y")
    figure
    hold on
end
for i = 1:num_iterations
    
    [~, ~, ~, ~, G]  = Intensity_to_F_to_PDF_calc(Q, R, I, f_sqrd, f_av_sqrd, I_compton, N_coh, N_inc, process);
    
    [G_slope, min_i] = Find_peak_minimum(R, G, guess_of_first_peak);
    
    
    if strcmp(process.MC_type_flag(2), "STD")
        std_current      = std(G(1:min_i) - G_slope);
    elseif strcmp(process.MC_type_flag(2), "INT")
        std_current      = trapz(R(1:min_i), abs(G(1:min_i) - G_slope));
    end
    
    % If better std, store the new parameters and decrease the search
    % interval
    if std_current < std_ref
        std_ref        = std_current;
        N_coh_ref      = N_coh;
        N_inc_ref      = N_inc;
        if strcmp(process.MC_type_flag(4), "span_VARY")
            span_interval  = span_interval/2;
        end
        if strcmp(process.plot_flag(1), "y")
            plot(R, G)
            plot(R(1:min_i), G_slope, '--k')
        end
    end
    
    N_coh = abs( span_interval(1)*(rand - 0.5) + N_coh_ref );
    
    if strcmp(process.MC_type_flag(3), "incoh_VARY")
        if or(strcmp(process.MC_type_flag(1), "spur_Q"), strcmp(process.MC_type_flag(1), "spur_Q2"))
            N_inc = span_interval(2)*(rand - 0.5)./10  + N_inc_ref ;
        else
            N_inc = abs(span_interval(2)*(rand - 0.5) + N_inc_ref );
            %N_inc = span_interval*(rand - 0.5) + N_inc_ref ;
        end
    end
    
    if mod(i, round(num_iterations/10)) == 0
        fprintf([num2str(100*i/num_iterations), ' percent done \n'])
    end
    
end
if strcmp(process.plot_flag(1), "y")
xlabel('R')
ylabel('G(R)')
hold off
pause(0.1)
end


end



