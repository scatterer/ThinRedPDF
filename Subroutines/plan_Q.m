clc
%close all;
clear;
format long g
% Code that determines the range in 2theta, step size, number of repeating
% cycles of a certain step size and measuring time for optimized PDF
% measurements. The color scheme in Fig 2, indicate the measuring times of
% that particular 2theta range


meas   = struct();
lambda = 0.7107488;

% Measuring time, tth range subdivided by length of "t"
t      = [542, 542, 776, 1383, 4500, 4500];
color  = ['c', 'c', 'm', 'b' ,  'r', 'r' ];

% Define Q range and largest acceptible 2theta step size 
Qmax        = 12;
Qmin        = 0.25;
dtth_accept = 0.02;

% Define starting Q step size and the final Q step size
dQ_min   = 0.04;
dQ_max   = 8*dQ_min;

% Variable dQ range. Create linear funtion dQ(Q) = dQ_slope*Q + dQ_off, that
% ranges between [dQ_min, dQ_max]. Invert and determine Q
dQ_slope = (dQ_max - dQ_min)/(Qmax - Qmin);
dQ_off   =  dQ_min - dQ_slope*Qmin;

% Initialize
Q(1) = Qmin;
i    = 2;
% While created Q vector is smaller than supplied Q max, determine slope
% and predictor for step size of index "i" and Q(i)
while Q(end) <= Qmax
    diff_dQ = 10;
    dQ_Q_est = dQ_slope*Q(i-1) + dQ_off;
    % While predictor deviate more than 1e-4 from linear function, create
    % new predictor
    while diff_dQ > 1e-4
        Q_est    = Q(i-1) + dQ_Q_est;
        dQ_Q_ref = dQ_Q_est;
        dQ_Q_est = dQ_slope*Q_est + dQ_off;
        diff_dQ  = abs(dQ_Q_ref - dQ_Q_est);
    end
    % Store final predictor Q(i) value
    Q = [Q, Q_est];
    i = i + 1;
end

% Determine the 2theta values, the final range of steps in Q and the
% corresponding step size in 2theta 
tth  = 2*asind(Q./4/pi*lambda);
dQ   = dQ_slope.*Q + dQ_off;
dtth = dQ./4/pi*lambda*180/pi*2./cosd(tth./2);

% Initialize the starting values for calculating the parameters for the
% measurement
meas.tth_start(1) = tth(1);
meas.tth_end(1)   = tth(1);
meas.dtth(1)      = dtth(1);
meas.cum_tth      = tth(1);
meas.cum_dtth     = dtth(1);
meas.t            = t(1);
meas.color        = color(1);

% Index
i   = 1;
c_1 = 1;
c_2 = 2;
c_t = 1;

% While the cumulative 2theta value is less than the final 2theta value,
% determine how many repeats of a certain step size is allowed with the set
% acceptible 2theta step size (dtth_accept)
while meas.cum_tth(end) < tth(end-1)
    dtth_diff = ( dtth(c_2) - dtth(c_1) );

    % While below the acceptible 2theta step size, increment in 2theta
    % until the criterion is breached
    while dtth_diff <= dtth_accept
        c_2       = c_2 + 1;
        dtth_diff = ( dtth(c_2) - dtth(c_1) );
    end

    % Count up the number of repeats of this particular step size, the
    % average step size, where the interval starts, where it ends and the
    % cumulative generated 2theta range and steps
    meas.counter(i)     = c_2 - c_1;
    meas.dtth(i)        = mean(dtth(c_1:c_2));
    meas.tth_start(i)   = meas.cum_tth(end);
    for j = 1:meas.counter(i)
        meas.cum_tth        = [meas.cum_tth, meas.cum_tth(end) + meas.dtth(i)];
        meas.cum_dtth       = [meas.cum_dtth, meas.dtth(i)];
    end
    meas.tth_end(i)     = meas.cum_tth(end);
    meas.tth_start(i+1) = meas.tth_end(i) + meas.dtth(i);

    % Determine if the 2theta range has gone passed a 1/length(t) portion
    % of the range. If so, change the measuring time accordingly
    if meas.tth_end(i) > c_t*tth(end)/length(t)
        c_t = c_t + 1;
    end
    meas.t(i)     = t(c_t);
    meas.color(i) = color(c_t);

    c_1 = c_2;
    i   = i + 1;
end

% With determined 2theta range, convert it to Q
meas.Q = 4*pi/lambda*sind(meas.cum_tth/2);



% figure
% plot(Q, tth,'k.');
% ylabel('$2\theta$', 'Interpreter','latex')
% xlabel('$Q$', 'Interpreter','latex')


figure
plot(Q, dQ, '.k')
ylabel('$\Delta Q$', 'Interpreter','latex')
xlabel('$Q$', 'Interpreter','latex')

figure
hold on
plot(tth, dtth, 'k.');
%plot(meas.cum_tth, meas.cum_dtth, 'r.');
for i = 1:length(meas.dtth)
    plot([meas.tth_start(i), meas.tth_end(i)],  meas.dtth(i).*[1,1], '-', 'Color', meas.color(i))
end
ylabel('$\Delta 2\theta$', 'Interpreter','latex')
xlabel('$2 \theta$', 'Interpreter','latex')

figure
hold on
dQ_est = dtth.*(2*pi/lambda)*(pi/180).*cosd(tth./2);
plot(4*pi/lambda*sind(0.5.*tth), dQ_est, 'k.');
%plot(meas.cum_tth, meas.cum_dtth, 'r.');
for i = 1:length(meas.dtth)
    dQ_steps =  meas.dtth(i).*(2*pi/lambda)*(pi/180).*cosd([meas.tth_start(i), meas.tth_end(i)]./2);
    plot(4*pi/lambda*sind(0.5.*[meas.tth_start(i), meas.tth_end(i)]),  dQ_steps, '-', 'Color', meas.color(i))
end
ylabel('$\Delta Q$', 'Interpreter','latex')
xlabel('$Q$', 'Interpreter','latex')

c     = 0;
cum_t = 0;
for i = 1:length(meas.dtth)
    c     = c + meas.counter(i);
    cum_t = cum_t + meas.t(i)*meas.counter(i);
    %disp(strcat(num2str(i), ", Repeat: ", num2str(meas.counter(i)), ", meas time: ", num2str(meas.t(i)), ", tth: ", num2str(round(meas.tth_start(i),3)), " - ", num2str(round(meas.tth_end(i),3)), ", step size: ", num2str(meas.dtth(i),4)))
    disp(strcat(num2str(i), ", Repeat: ", num2str(meas.counter(i)), ", meas time: ", num2str(meas.t(i)./542), " * 542 s, tth: ", num2str(round(meas.tth_start(i),3)), " - ", num2str(round(meas.tth_end(i),3)), ", step size: ", num2str(meas.dtth(i),4)))
end
disp(' ')
disp(strcat("Total number of steps: ", num2str(c), ", Total measuring time: ", num2str(cum_t/3600), " h" ))
