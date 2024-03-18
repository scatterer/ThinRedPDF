function [a,da] = linfit(x,y,e)
%----------------------------------------------------------
% LINFIT    Straight line fitting of data width errors
%           in the y-coordiante.
%
% [a,da] = linfit(x,y,e)
%
% x  - a vector of x-data values
% y  - a vector of y-data values
% e  - a vector of errors associated with y-data values
%      (optional parameter)
%
% a  - fitted line parameters: y(x) = a(1)*x + a(2)
% da - variances of fitted line parameters
%
% example: a linear regression
%
%              x = [0:10];
%              y = 2.0*x + 3;
%              r = rand(size(x))-0.5; % random deviations
%              y = y + r;
%              e = std(r)*ones(size(x)); % errors
%              a = linfit(x,y,e);
%              errorbar(x,y,e,'k*')
%              hold on, plot(x,polyval(a,x),'b')
%
% References:
%  [1] W.H.Press, S.A.Teukolsky, W.T.Vetterling,
%      B.P.Flannery, NUMERICAL RECIPES IN C
%
% See also POLYVAL, POLYFIT
%
% version 1.0, 14.3.2005, (c) Zdenek Matej
%----------------------------------------------------------

if ~exist('e','var') | isempty(e) e = ones(size(x)); end
x = x(:); y = y(:); e = e(:);
S = 1./(e.*e); Sx = x'*S; Sy = y'*S; S = sum(S);
t = 1./e.*(x-Sx/S); Stt = t'*t;
a(1)  = 1/Stt*t'*(y./e);
a(2)  = (Sy-Sx*a(1))/S;
da(1) = sqrt(1/Stt);
da(2) = sqrt(1/S*(1+Sx^2/S/Stt));
return;
