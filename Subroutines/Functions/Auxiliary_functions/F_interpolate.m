function [Q_inter, F_inter] = F_interpolate(Q, F)


coeff    = polyfit([0; Q(1:5)], [0; Q(1:5).*F(1:5)], 1);
dQ       = Q(2) - Q(1);
Q_init   = transpose([1e-20:dQ:Q(1)]);
F_Q_init = polyval(coeff, Q_init);

Q_inter = [Q_init; Q(2:end)];
F_inter = [F_Q_init./Q_init; F(2:end)];
F_inter(1) = F_inter(2);
end




