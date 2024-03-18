close all;
clear all;

thedir  = './data/batch_static/';
thedir  = './data/batch_static_pc/';
thedir  = './data/batch_static_pc_omega018/';



files = dir([thedir,'*.xrdml']);

file_open = './data/open_detector_1.xrdml';
file_open_pc = './data/open_detector_parallel_collimator_1.xrdml';

d      = XRDMLread(file_open);
tth_open = d.x;
I_open   = d.data;

d      = XRDMLread(file_open_pc);

tth_open_pc = d.x;
I_open_pc   = d.data;


R = 240;
omega = 0.375;

D = ((1:1:127)-64)*8.89/127;
tth_expected = 1.94 + 2*atand(D/2/R);

F = 20;

for i = 1:length(files)
  d          = XRDMLread([thedir,files(i).name]);
  tth{i}     = d.x;
  
  pixel_positions = 2*R*tand( (d.x-median(d.x))/2 );
  
  beta(i) = median(d.x) - omega;
  
  % Size of diffracted beam if only accepting parallel rays (using a
  % perfect collimator)
  
  Diffracted_size(i) = F*sind(beta(i));
  
  Detector_footprint(i) = pixel_positions(end)*2/sind(beta(i));
  
  % Largest possible footprint that the detector can sustain over the
  % entire 2theta range is 8.89 mm. It is OK that this value is lower than
  % the illuminated beam, since it varies as sin(beta) we get a constant
  % factor, instead of a nasty sind(beta).
  
  Detector_size_constant_footprint(i) = pixel_positions(end)*2*sind(beta(i));
  
  [pos_diffracted_size,Id_diffracted_size] = find(pixel_positions >= -Diffracted_size(i)/2 & pixel_positions <= Diffracted_size(i)/2);
  
  [pos_dscf,Id_dscf] = find(pixel_positions >= -Detector_size_constant_footprint(i)/2 & pixel_positions <= Detector_size_constant_footprint(i)/2);
  
  
  effective_detector_positions_dscf = pixel_positions(Id_dscf);
  effective_detector_opening_dscf(i) = abs(effective_detector_positions_dscf(1)) + effective_detector_positions_dscf(end);

  effective_detector_positions_ds = pixel_positions(Id_diffracted_size);
  effective_detector_opening_ds(i) = abs(effective_detector_positions_ds(1)) + effective_detector_positions_ds(end);

  I{i}         = d.data(Id_diffracted_size);
  Itot_diffracted_size(i)      = sum(d.data(Id_diffracted_size));
  Itot_dscf(i)  = sum(d.data(Id_dscf));
  tth_tot(i)   = median(d.x);
  
end

I_open_ip = interp1(tth_open,I_open,tth_tot);
I_open_pc_ip = interp1(tth_open_pc,I_open_pc,tth_tot);


fig = figure('Color','w');
plot(tth_tot,Itot_diffracted_size/max(Itot_diffracted_size),'ko');
hold on;
plot(tth_tot,Itot_dscf/max(Itot_dscf),'ks');
plot(tth_tot,(I_open_ip*1.04)/max(I_open_ip*1.04),'kx');
plot(tth_tot,I_open_pc_ip/max(I_open_pc_ip),'kd');

xlabel('2\theta [degrees]');
ylabel('Intensity [cps]');
set(gca,'FontSize',14);
axis([0 90 0 1.3]);
legend('Matching diffracted beam size','constant detector foot print','open detector','pc');

fig = figure('Color','w');
plot(tth_tot,Itot_diffracted_size,'ko');
hold on;
plot(tth_tot,Itot_dscf,'ks');
plot(tth_tot,I_open_ip*1.04,'kx');
plot(tth_tot,I_open_pc_ip,'kd');

xlabel('2\theta [degrees]');
ylabel('Intensity [cps]');
set(gca,'FontSize',14);
axis([0 90 0 Inf]);
legend('Matching diffracted beam size','constant detector foot print','open detector','pc');


figure('Color','w');
plot(beta,Diffracted_size,'kx')
hold on;
plot(beta,effective_detector_opening_ds,'ko')
plot(beta,effective_detector_opening_dscf,'ks')
legend('Diffracted size','Detector opening DS','Detector opening CF');
xlabel('\beta [degrees]');
ylabel('Width [mm]');
%plot(beta,Detector_footprint,'bo');

[tth_sorted,ID] = sort(tth_tot,'ascend');
Itot_dscf_sorted = Itot_dscf(ID);

fig = figure('Color','w');
hold on;
plot(tth_sorted,Itot_dscf_sorted,'-ko');

xlabel('2\theta [degrees]');
ylabel('Intensity [cps]');
set(gca,'FontSize',14);
axis([0 90 0 Inf]);
legend('constant detector foot print with collimator');

