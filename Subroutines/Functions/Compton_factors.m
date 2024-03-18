function [I_compton] = Compton_factors(Q, CC, composition, Z, wavelength)

form_factor = zeros(length(composition), length(Q));
f           = zeros(length(Q), 1);

for i = 1:length(composition)
    for j = 1:5
        gaussian_basis   = CC.a(i,j).*exp( -CC.b(i,j).*(Q/(4.*pi)).^2);
        form_factor(i,:) = form_factor(i,:) + transpose(gaussian_basis);
    end
end

for i=1:length(composition)
    f = f + transpose( composition(i).*form_factor(i,:));
end

mean_Z    = sum(composition.*Z);
I_compton = mean_Z.*ones(length(Q),1) - f;


%Icomptom from Thijsse
n = length(Z);
i = 0;

h      = 6.62607004e-34;
m      = 9.10938356e-31;
c0     = 299792458;

lambda = wavelength*1e-10;
theta  = asind(Q.*1e10/4/pi*lambda);
lambdap = lambda + 2*h/(m*c0).*sind(theta).^2;

c = composition;
a = 2.6917.*Z.^(-1) + 1.245;
b = 1.1870.*Z.^(-1) + 0.1075 + 0.00436.*Z - (0.01543.*Z).^2 + (0.01422.*Z).^3;

for j = 1:n
    i = i + c(j)*Z(j)*(b(j).*Q).^a(j)./(1+(b(j).*Q).^a(j));
end

%I_compton =  (lambda./lambdap).^2.*i;

I_compton = (lambda./lambdap).^2.*I_compton;
%I_compton = I_compton;
%(lambda./lambdap).^2
% figure
% plot(Q,I_compton_alt,'-o');
% hold on;
% plot(Q,I_compton,'-x');
% plot(Q,I_compton_alt./(lambda./lambdap).^2,'--x');
% legend('Thijsse', 'Johan','Thijsse no Breit')

end
