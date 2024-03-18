%------------------------------------------------------%
% example3.m   An exaple how XRDML files can be        %
%              loaded using the XRDMLread function.    %
%                                                      %   
% Sin^2(Psi) method applied to multiple scans          %
% of an anatase sample.                                %   
%                                                      %
% version 1.1, 22.11.2010                              %
% (c) Zdenek Matej, Milan Dopita                       %
%------------------------------------------------------%

d = XRDMLread('TiO2_STRESS.xrdml');

figure

hold on
col = hsv(length(d.Psi));
for k=1:length(d.Psi)
    plot3( d.Theta2(k,:) , repmat(d.Psi(k),1,size(d.data,2)) , ...
           round(d.data(k,:)*d.time), '*' , 'Color' , col(k,:) )
end
view(45,20)

xlabel( [d.xlabel ' (' d.xunit ')'] )
ylabel( [d.ylabel ' (' d.yunit ')'] )
zlabel( 'Intensity (counts)' )

title( d.filename , 'Interpreter','none','Interpreter','none','FontName','Helvetica')

% set global wavelength data for the pseudoVoigt function
global WAVELENGTHS
switch d.kType
case {'K-Alpha 1'}
    WAVELENGTHS = [1.0 0.0];
case {'K-Alpha'}
    WAVELENGTHS = [           1.0 (d.kAlpha1-d.Lambda)/d.Lambda ;
                    d.kAlphaRatio (d.kAlpha2-d.Lambda)/d.Lambda ];
otherwise
    warning('usedWavelength type not supported (using K-Alpha 1)')
    WAVELENGTHS = [1.0 0.0];
end

% add diffraction lines (anatase - database)
a0(1,:) = [56 36.95 0.2 0.5]; % (103)
a0(2,:) = [88 37.79 0.2 0.5]; % (004)
a0(3,:) = [36 38.57 0.2 0.5]; % (112)

% set refined lines parameters
Linda(1,:) = [1 1 1 0]; % (103)
Linda(2,:) = [1 1 1 1]; % (004)
Linda(3,:) = [1 1 1 0]; % (112)

results(1).hkl = [1 0 3];
results(2).hkl = [0 0 4];
results(3).hkl = [1 1 2];

% add linear background
b0 = [0.0 0.0];

% fit all scans
for k=1:length(d.Psi)
    x = d.Theta2(k,:);
    y = round(d.data(k,:)*d.time);
    w = 1./y;
    [a,b,da,db] = pseudoVoigtFit(x,y,w,a0,b0,[],Linda);
    yc = sum(pseudoVoigt(a,x),1) + polyval(b,x);
    for m=1:length(results)
        results(m).T2(k,1)  =  a(m,2);
        results(m).dT2(k,1) = da(m,2);
    end
    plot3( x , repmat(d.Psi(k),1,size(d.data,2)) , yc , 'k' )
end

% calc lattice planes distaces and their deviations
for m=1:length(results)
    th = results(m).T2*pi/360;
    results(m).d  = d.Lambda./2./sin(th);
    results(m).dd = cos(th)/d.Lambda.*results(m).d.^2.* ...
               results(m).dT2*pi/180;
end

% create sin^2(Psi) plots
x = sin(d.Psi*pi/180).^2;

for m=1:length(results)
    
    figure
    
    y = results(m).d; e = results(m).dd;
    [results(m).p,results(m).dp] = linfit(x,y,e);
    if ~exist('OCTAVE_VERSION'), %#ok<EXIST>
        errorbar( x , y , e , 'k*' ), hold on
    else
        plot( x , y , 'k*' ), hold on,
    end
    plot( x , results(m).p(1)*x+results(m).p(2) , 'k' )
    xlim( [-0.03, 0.9] )
    xlabel( 'sin^2(Psi)' ), ylabel ( 'd (A)' )
    str = sprintf('ANATASE (%d%d%d)', results(m).hkl);
    title( str )
end

% print results
for m=1:length(results)
    str = sprintf('   ANATASE DIFFRACTION LINE (%d%d%d)\n', ...
                  results(m).hkl);
    str = [ str sprintf('%14s%26s\n','Psi (deg)','2Theta (deg)') ];
    str = [ str sprintf('%14.4f%15.4f +/- %6.4f\n', ...
                 [d.Psi results(m).T2 results(m).dT2]' ) ];
    str = [ str sprintf('   D(HKL) VS. SIN^2(PSI)\n') ];
    str = [ str sprintf('%31s%23s\n','slope','intercept (A)') ];
    str = [ str sprintf('%15.4e +/- %11.4e %10.5f +/- %7.5f\n', ...
                 [results(m).p; results(m).dp]) ];
    disp(str)
end
