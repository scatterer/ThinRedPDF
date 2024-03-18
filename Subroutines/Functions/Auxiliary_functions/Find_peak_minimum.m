
function [G_slope, min_i] = Find_peak_minimum(R, G, R_guess)


%[~, max_guess_i] = min(abs(R - R_guess));

[~, lower_bound] = min(abs(R - (R_guess - 0.5)));
[~, upper_bound] = min(abs(R - (R_guess + 0.5)));

%lower_bound = max_guess_i - 10;
%upper_bound = max_guess_i + 10;
[~, max_i]  = max(G(lower_bound: upper_bound));
max_i       = max_i + lower_bound;



for k = 1:max_i
    delta_G = (G(max_i-k) - G(max_i-k-1) );
    if delta_G < 0
        min_i = max_i - k;
        break
    end
end
G_slope = transpose((G(min_i)/R(min_i)).*R(1:min_i));

end