function [x,y,t,e] = process_pdf(thedir,omega,apply_correction)
R = 240;
F = 20;

files = dir([thedir,'*.xrdml']);

[sort_name,ind]=natsortfiles({files.name});


for i = 1:length(files)
    if mod(i,75) == 0
        disp([num2str(100.*i/length(files)),  ' Percent done'])
    end
    %d          = XRDMLread([thedir,files(i).name]);
    d          = XRDMLread(strcat(thedir, sort_name{i}));

    %strcat(thedir, sort_name{i})


    if apply_correction
        %tth{i}     = d.x;

        pixel_positions = 2*R*tand( abs(d.x - median(d.x))./2 );

        %pixel_positions = 0.1*tand( abs(d.x - median(d.x))./2 );


        %beta(i) = median(d.x) - omega;
        beta = d.x + omega;

        % Size of diffracted beam if only accepting parallel rays (using a
        % perfect collimator)

        %Diffracted_size(i) = F*sind(beta(i));

        %Detector_footprint(i) = 8.89/sind(beta(i));

        % Largest possible footprint that the detector can sustain over the
        % entire 2theta range is 8.89 mm. It is OK that this value is lower than
        % the illuminated beam, since it varies as sin(beta) we get a constant
        % factor, instead of a nasty sind(beta).

        %8.89*sind(beta(i));
        Detector_size_constant_footprint = 8.89*sind(beta);


        %[pos_diffracted_size,Id_diffracted_size] = find(pixel_positions >= -Diffracted_size(i)/2 & pixel_positions <= Diffracted_size(i)/2);

        %[pos_dscf,Id_dscf] = find(pixel_positions > -Detector_size_constant_footprint/2 & pixel_positions < Detector_size_constant_footprint/2);
        [~,Id_dscf] = find( pixel_positions <= Detector_size_constant_footprint./2);


        %effective_detector_positions_dscf = pixel_positions(Id_dscf);
        %effective_detector_opening_dscf(i) = abs(effective_detector_positions_dscf(1)) + effective_detector_positions_dscf(end) + 0.07;

        %effective_detector_positions_ds = pixel_positions(Id_diffracted_size);
        %effective_detector_opening_ds(i) = abs(effective_detector_positions_ds(1)) + effective_detector_positions_ds(end) + 0.07;

        %I{i}         = d.data(Id_diffracted_size);
        %Itot_diffracted_size(i)      = sum(d.data(Id_diffracted_size));
        Itot_dscf(i) = sum(d.data(Id_dscf));
        %Itot_dscf(i) = sum(d.data);
        count(i)     = sum(d.data(Id_dscf).*d.time);


        I_matrix{i} = d.data(:);
        len_pixel(i) = length(d.data(Id_dscf));


    else
        Itot_dscf(i)  = sum(d.data);
    end
    uItot_dscf(i) = sqrt( Itot_dscf(i)./d.time ) ;
    time_dscf(i)  = d.time;
    tth_tot(i)    = median(d.x);
    t(i)          = length(d.data);
end
[tth_sorted,ID]   = sort(tth_tot,'ascend');
Itot_dscf_sorted  = Itot_dscf(ID);
uItot_dscf_sorted = uItot_dscf(ID);
meas_time         = time_dscf(ID);
len_pixel         = len_pixel(ID);
count             = count(ID);

for i = 1:length(I_matrix)
    I_matrix_sort{i}     = I_matrix{ID(i)};
end
x = tth_sorted(:);
y = Itot_dscf_sorted(:);
e = uItot_dscf_sorted(:);
t = meas_time(:);


dev = len_pixel'./(length(d.x).*sind(x));
%
% figure
% hold on
% plot(x, len_pixel,'ro')
% plot(x, length(d.x).*sind(x), 'ok')


% figure
% hold on
% step = 0;
% for i = 1:length(I_matrix)
%     step = step + 0.05;
% plot(I_matrix_sort{i}./max(I_matrix_sort{i}) + step)
% end
figure
hold on
pcolor([1:127], x, cell2mat(I_matrix_sort)')
plot(len_pixel./2 + 127/2, x, 'k')
plot(-len_pixel./2 + 127/2, x, 'k')
shading flat
set(gca,'ColorScale','log')


% Imatrix = struct();
% Imatrix.x = x;
% Imatrix.len = [1:127];
% Imatrix.len_pixel = len_pixel;
% Imatrix.matrix = cell2mat(I_matrix_sort)';
% 
% save('I_matrix.mat', 'Imatrix')
% return

I_matrix_sort = cell2mat(I_matrix_sort)';
c = 1;
i = 1;
while i < length(y) - 1
    if x(i+1) == x(i)
        y_dummy = 0;
        e_dummy = 0;
        t_dummy = 0;
        counter = 0;
        I_matrix_dummy = zeros(1, 127);

        while x(i+1) == x(i)
            if i + 1 == length(y)
                break
            end
            counter = counter + 1;
            y_dummy(counter) = y(i);
            e_dummy(counter) = e(i).^2;
            t_dummy(counter) = t(i);
            I_matrix_dummy   = I_matrix_dummy + I_matrix_sort(i,:);
            
            i = i + 1;
        end
        counter = counter + 1;
        y_dummy(counter) = y(i);
        e_dummy(counter) = e(i).^2;
        t_dummy(counter) = t(i);
        I_matrix_dummy   = I_matrix_dummy + I_matrix_sort(i,:);
        %         y(i-counter:i)
        %         y_dummy
        %         x(i-counter:i)
        y_repeat(c) = mean(y_dummy);

        e_repeat(c) = sqrt(sum(e_dummy))./length(e_dummy);
        t_repeat(c) = sum(t_dummy);
        I_matrix_repeat(c,:) = I_matrix_dummy;
    else
        y_repeat(c) = y(i);
        e_repeat(c) = e(i);
        t_repeat(c) = t(i);
        I_matrix_repeat(c,:) = I_matrix_sort(i,:);
    end
    x_repeat(c) = x(i);
    len_pixel_repeat(c) = len_pixel(i);
    c = c + 1;
    i = i + 1;
end


figure
hold on
pcolor([1:127], x_repeat, I_matrix_repeat)
plot(len_pixel_repeat./2 + 127/2, x_repeat, 'k')
plot(-len_pixel_repeat./2 + 127/2, x_repeat, 'k')
shading flat
set(gca,'ColorScale','log')

figure
hold on
for i = 1:length(x_repeat)
    plot([1:127], I_matrix_repeat(i,:) + max(I_matrix_repeat(i,:)) .*i, 'k-')
end

figure
hold on
plot(x_repeat, y_repeat, 'ok')
plot(x, y, '.')

x = x_repeat';
y = y_repeat';
e = e_repeat';
t = t_repeat';

% figure
% plot(x, t)
% figure
% hold on
% plot(x, e, 'rx')
% plot(x, sqrt(y./t), 'ko')

% max_v
% set(gca,'colorscale','log')
% ylabel('tth [deg]')
% xlabel('Pixel [mm]')
% pause()

end