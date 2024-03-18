fig = figure();
set(gcf,'Color','w')
set(gca, 'YTick', []);
set(gca, 'XTick', []);
set(gca,'XColor', 'none','YColor','none')
warning('off')

[ha, pos] = tight_subplot(3,2,[.11 .08],[.08 .08],[.045 .01]);

ax2 = axes('position', pos{1});
ax2.YColor = 'none';
ax2.Color = 'none';

Q_vec = zeros(1,length([0:10:data.tth(end)+15]));
for i = 1:length([0:10:data.tth(end)+15])
    Q_vec(i) = 4*pi/instr.wavelength.*sind((i-1)*10./2);
end

ax2.XTick = Q_vec;
xlim(ax2, [0, Q_vec(end)])
xtickformat('%.1f')
set(gca, 'TickLength', [0.01, 0.02])
ax2.XAxisLocation = 'top';
xlabel(ax2, '$Q$ [\AA$^{-1}$]', 'Interpreter','latex')
set(gca,'FontSize',13)

axes(ha(1))
hold on
plot(data.tth, data.I_scan, 'k.-',  'DisplayName','Scan')
errorbar(data.tth, data.I_scan, data.I_error, 'k.', 'HandleVisibility','off')

if strcmp(process.bkgr_flag, 'y')
    plot(theta_bkgr, I_bkgr, 'r', 'DisplayName','Background')
end
plot(data.tth, data.I_corr_orig, 'r.-',  'DisplayName','Corrected scan')
errorbar(data.tth, data.I_corr_orig, data.I_corr_error, 'r.', 'HandleVisibility','off')

hold off
set(gca, 'TickLength', [0.02, 0.02])
set(gca, 'XTick', 0:10:data.tth(end)+15)
xlabel('$2\theta$ [deg]', 'Interpreter','latex')
ylabel('$I$ [cps]', 'Interpreter','latex')
set(gca,'FontSize',13)
%set(gca, 'Box', 'on')
ytickformat('%.0f')
xtickformat('%.0f')
legend

% ---------------------------------------------------------------------------------------------------------------------------------------

ax2 = axes('position', pos{2});
ax2.YColor = 'none';
ax2.Color = 'none';

tth_vec = zeros(1,length([0:1:data.Q(end)+1]));
for i = 1:length([0:1:data.Q(end)+1])
    tth_vec(i) = 2.*asind((i-1).*instr.wavelength./(4*pi));
end

ax2.XTick = tth_vec;
xlim(ax2, [0, tth_vec(end)])
xtickformat('%.0f')
set(gca, 'TickLength', [0.01, 0.02])
ax2.XAxisLocation = 'top';
xlabel(ax2, '$2\theta$ [deg]', 'Interpreter','latex')
set(gca,'FontSize',13)

axes(ha(2))
plot(data.Q, data.I_norm, '-k'    )
hold on
errorbar(data.Q, data.I_norm, data.I_norm_error, 'k.', 'HandleVisibility','off')
plot(data.Q, data.f_sqrd, '--b'   )
plot(data.Q, data.f_av_sqrd, '--', 'color', '#336600')
plot(data.Q, data.N_inc(1).*data.I_compton.*data.corr.geo(1:length(data.I_compton)), '--r')
if strcmp(sample.sub.flag, "y")
    plot(data.Q, data.N_inc(2).*data.sub.I_compton.*data.corr.geo(1:length(data.I_compton)), '-.r')
end
hold off
xlim([0,data.Q(end)+1])
set(gca, 'TickLength', [0.01, 0.02])
set(gca, 'XTick', -0.1:1:data.Q(end)+1)
set(gca, 'YAxisLocation', 'left')
xlabel('$Q$ [\AA$^{-1}$]', 'Interpreter','latex')
ylabel('I normalized [cps]', 'Interpreter','latex')
set(gca,'FontSize',13)
ytickformat('%.0f')
xtickformat('%.0f')
set(gca, 'Box', 'off')
legend(["Normalized I";"$<f^2>$"; "$<f>^2$"; "Compton I"; "Sub compton I"],  'location', 'northeast', 'interpreter', 'latex')



% ---------------------------------------------------------------------------------------------------------------------------------------

SQ = load('g_G_F_I_V33Zr67_Hx_no_window');


ax2 = axes('position', pos{3});
ax2.YColor = 'none';
ax2.Color = 'none';
ax2.XTick = tth_vec;
xlim(ax2, [0, tth_vec(end)])
xtickformat('%.0f')
set(gca, 'TickLength', [0.01, 0.02])
ax2.XAxisLocation = 'top';
xlabel(ax2, '', 'Interpreter','latex')
set(gca,'FontSize',13)

axes(ha(3))
hold on
errorbar(data.Q, data.Q.*data.F_orig, data.F_error, 'k.', 'HandleVisibility','off')
plot(data.Q, data.Q.*data.F_orig, 'k.')
plot(data.Q_inter, data.Q_inter.*data.F, 'r-')
plot(data.Q_inter, zeros(length(data.Q_inter)), 'k', 'HandleVisibility','off')
plot(SQ.Q, SQ.H00.F_weight, 'b-')
hold off
xlim([0,data.Q_inter(end)+1])
ylim([min(data.Q.*data.F_orig) - 1, max(data.Q.*data.F_orig) + 1])
set(gca, 'TickLength', [0.02, 0.06])
set(gca, 'XTick', 0:1:data.Q(end)+1)
xlabel('$Q$ [\AA$^{-1}$]', 'Interpreter','latex')
ylabel('$Q\cdot F$', 'Interpreter','latex')
set(gca,'FontSize',13)
%set(gca, 'Box', 'on')
ytickformat('%.0f')
xtickformat('%.0f')
legend(["F", "SQ"], 'interpreter', 'latex', 'Location', 'northeast')


% ---------------------------------------------------------------------------------------------------------------------------------------

axes(ha(4))
%plot(data.R, data.G,'-k.')
hold on
plot(data.R, data.G_L,'-r')
%plot(data.R, data.G' + data.R,'-k.')
plot(data.R, zeros(length(data.R),1), '--k')
%plot(data.R, data.R'.*ones(length(data.R),1), '--k')
plot(SQ.H00.R, SQ.H00.G_weight, 'b.')
ylabel('G(R)')
xlabel('R')
hold off
xlim([0,data.R(end)])
set(gca, 'TickLength', [0.02, 0.06])
set(gca, 'XTick', 0:1:data.R(end))
xlabel('$R$ [\AA]', 'Interpreter','latex')
ylabel('$G(R)$', 'Interpreter','latex')
set(gca,'FontSize',13)
set(gca, 'Box', 'on')
ytickformat('%.0f')
xtickformat('%.0f')

% ---------------------------------------------------------------------------------------------------------------------------------------

axes(ha(5))
hold on
plot(SQ.H00.R, SQ.H00.g_weight, 'b.')
%plot(R_warren, g_warren, 'b.')
%plot(data.R, data.g,'-k')
plot(data.R, data.g_L,'-r')
plot(data.R, ones(length(data.R),1), '--k')
ylabel('g(R)')
xlabel('R')
hold off

[~, min_ind] =min(abs(data.R -1.5));
xlim([0,data.R(end)])
ylim([-1, max(data.g(min_ind:end)) + 0.5])
set(gca, 'TickLength', [0.02, 0.06])
set(gca, 'XTick', 0:1:data.R(end))
xlabel('$R$ [\AA]', 'Interpreter','latex')
ylabel('$g(R)$', 'Interpreter','latex')
set(gca,'FontSize',13)
set(gca, 'Box', 'on')
ytickformat('%.0f')
xtickformat('%.0f')

% ---------------------------------------------------------------------------------------------------------------------------------------

axes(ha(6))
plot(data.tth(1:length(data.corr.geo)), data.corr.geo, '--', 'color', '#336600', 'DisplayName','E-window correction')
hold on
plot(data.tth, data.corr.pol, '--', 'color', '#cc9900', 'DisplayName','Polarization correction')
plot(data.tth, data.corr.abs, '--b', 'DisplayName','Absorption correction')
plot(data.tth, data.corr.tot, 'k', 'DisplayName','Total correction')

plot(data.tth(1:length(data.I_compton)), 0.5.*data.N_inc(1).*data.I_compton.*data.corr.geo(1:length(data.I_compton))./max(data.N_inc(1).*data.I_compton.*data.corr.geo(1:length(data.I_compton))), '--r', 'DisplayName', 'Compton profile')
% if strcmp(sample.sub.flag, "y")
%     plot(data.tth, 0.5.*data.N_inc(2).*data.sub.I_compton.*data.corr.geo./max(data.N_inc(2).*data.sub.I_compton.*data.corr.geo), '-.r', 'HandleVisibility','off')
% end

hold off
ylim([0.0,1.05])
set(gca, 'TickLength', [0.02, 0.06])
set(gca, 'XTick', 0:15:data.tth(end)+15)
xlabel('$2\theta$ [deg]', 'Interpreter','latex')
ylabel('Corrections', 'Interpreter','latex')
set(gca,'FontSize',13)
ytickformat('%.1f')
xtickformat('%.0f')
legend('Position', [0.48, 0.28, 0.1, 0.1],'FontSize',12)


function [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   subplot
%        pos    positions of the axes objects
%
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
% Pekka Kumpulainen 21.5.2012   @tut.fi
% Tampere University of Technology / Automation Science and Engineering
if nargin<3; gap = .02; end
if nargin<4 || isempty(marg_h); marg_h = .05; end
if nargin<5; marg_w = .05; end
if numel(gap)==1;
    gap = [gap gap];
end
if numel(marg_w)==1;
    marg_w = [marg_w marg_w];
end
if numel(marg_h)==1;
    marg_h = [marg_h marg_h];
end
axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh;
axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;
py = 1-marg_h(2)-axh;
% ha = zeros(Nh*Nw,1);
ii = 0;
for ih = 1:Nh
    px = marg_w(1);
    
    for ix = 1:Nw
        ii = ii+1;
        ha(ii) = axes('Units','normalized', ...
            'Position',[px py axw axh], ...
            'XTickLabel','', ...
            'YTickLabel','');
        px = px+axw+gap(2);
    end
    py = py-axh-gap(1);
end
if nargout > 1
    pos = get(ha,'Position');
end
ha = ha(:);
end
