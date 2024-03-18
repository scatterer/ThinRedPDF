
function Coordination_num = Calc_coordination_number(R, g, num_density, R_guess)


[~, max_guess_i] = min(abs(R - R_guess));

g_temp = g - 1;

lower_bound = max_guess_i - 10;
upper_bound = max_guess_i + 10;
[~, max_i]  = max(g_temp(lower_bound: upper_bound));
max_i       = max_i + lower_bound;



for k = 1:max_i
    if g_temp(max_i - k) < 0
        lower_bound = max_i - k;
        break
    end
end

for k = 1:max_i
    if g_temp(max_i + k) < 0
        upper_bound = max_i + k;
        break
    end
end

R_peak = transpose( R(lower_bound:upper_bound) );
g_peak = g(lower_bound:upper_bound);


Coordination_num = 4.*pi.*num_density.*trapz(R_peak, R_peak.^2.*(g_peak - 1) );


end