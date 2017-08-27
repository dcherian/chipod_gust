function [] = DebugPlots(hfig, t0, t1, chi, name, ww)

    if ~exist('ww', 'var'), ww = 1; end
    if isempty(hfig), hfig = gcf(); end

    i0 = find_approx(chi.time, t0, 1);
    i1 = find_approx(chi.time, t1, 1);
    tind = i0:i1;

    figure(hfig);

    time = moving_average(chi.time(tind), ww, ww);

    try
        figure(hfig);
    catch ME
        CreateFigure;
    end
    ax(1) = subplot(511);
    semilogy(time, moving_average(chi.chi(tind), ww, ww), 'displayname', name)
    ylabel('\chi')
    Common()

    ax(2) = subplot(512);
    plot(time, moving_average(chi.dTdz(tind), ww, ww), 'displayname', name)
    hold on;
    plot(xlim, [0, 0], 'k--');
    ylabel('dT/dz')
    Common()
    % symlog(gca, 'y', 5e-3);

    ax(3) = subplot(513);
    try
        semilogy(time, moving_average(chi.eps(tind), ww, ww), 'displayname', name)
    catch ME
        semilogy(time, moving_average(chi.eps1(tind), ww, ww), 'displayname', name)
    end
    ylabel('\epsilon')
    ylim([10.^[-10, -3]])
    Common()

    ax(4) = subplot(514);
    try
        semilogy(time, moving_average(chi.Kt(tind), ww, ww), ...
                 'displayname', name)
    catch ME
        semilogy(time, moving_average(chi.Kt1(tind), ww, ww), ...
                 'displayname', name)
    end
    ylabel('K_t')
    ylim([10.^[-7, 0]])
    Common()

    ax(5) = subplot(515);
    try
        plot(time, moving_average(chi.Jq(tind), ww, ww), 'displayname', name)
    catch ME
        plot(time, moving_average(-chi.Jq1(tind), ww, ww), 'displayname', name)
    end
    ylabel('J_q^t')
    Common()
    legend('-DynamicLegend');

    linkaxes(ax, 'x')
    xlim([t0, t1])
    datetick('x', 'mmm-dd HH:MM', 'keeplimits')
end

function Common()
    hold on
    datetick('x', 'mm-dd HH:MM', 'keeplimits')
end