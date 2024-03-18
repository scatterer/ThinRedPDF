
function G = PDF_calculation(Q, R, F)

G = zeros(length(R), 1);
for i = 1:length(R)
    integrand = Q.*F.*sin(Q.*R(i));
    G(i) = trapz(Q, integrand);
end

G = 2.*G./pi;



