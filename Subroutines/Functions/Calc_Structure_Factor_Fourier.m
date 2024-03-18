function [Q, F, S] = Calc_Structure_Factor_Fourier(r, G, Q, trunc_flag, guess_of_first_peak)

%Q = transpose(linspace(0, 50,1000));
F = zeros(length(Q),1);


%Landczos = sin(pi.*r./max(r))./(pi.*r./max(r));
%G        = G.*Landczos;


if trunc_flag == 1
    [G_slope, min_i] = Find_peak_minimum(r, G, guess_of_first_peak);

    r = r(1:min_i);
    G = G(1:min_i) - G_slope';
    %G = G(1:min_i);
   
    %Landczos = sin(pi.*r./max(r))./(pi.*r./max(r));
    %G        = G.*Landczos;
end
for i = 1:length(Q)
    F(i) = trapz(r, G.*sin(Q(i).*r));
end

F = F./Q;
S = F + 1;
end