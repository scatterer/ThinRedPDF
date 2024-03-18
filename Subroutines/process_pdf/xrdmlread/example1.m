%------------------------------------------------------%
% example1.m   An exaple how XRDML files can be        %
%              loaded using the XRDMLread function.    %
%                                                      %  
% Demonstration of reading a single Omega/2Theta scan  %
%                                                      %
% version 1.1, 22.11.2010                              %
% (c) Zdenek Matej, Milan Dopita                       %
%------------------------------------------------------%

% load XRDML data by a simple single command      
d = XRDMLread('BST_O2T.xrdml') %#ok<NOPTS>

% plot data
figure

semilogy( d.x , round(d.data*d.time) )

box on

xlabel( [d.xlabel '(' d.xunit ')'] )
ylabel( 'Intensity (counts)' )

title( d.filename , 'Interpreter','none','FontName','Helvetica')

