function rho0 = sum_rule(F,Q)
   % sum(Qi*F(Qi)*dQ) = -2*pi^2*rho0
   % number of atoms/A^3
   % Johans F does not include Q factor, so we need to multiply by Q^2
   % int([A^-2]*[A-2]) = A^-3
   rho0 = -1/(2*pi^2)*trapz(Q,Q.^2.*F);
end