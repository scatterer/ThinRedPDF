
function F = F_fourier_calculation(Q, R, G)

F = zeros(length(Q), 1);
size(G)
for i = 1:length(Q)
    integrand = R.*G.*sin(Q(i).*R);
    F(i) = trapz(R, integrand)./(4*pi);
end



