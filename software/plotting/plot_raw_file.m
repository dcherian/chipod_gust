function plot_raw_file(fid)

    here    =   pwd;                % mfiles folder
    basedir =   here(1:(end-6));    % substract the mfile folder
    head = chi_get_calibration_coefs(basedir);

    [data, ~] = raw_load_chipod(fid);
    cal = chi_calibrate_chipod(data, head)

    rawdate = data.datenum(1:2:end);
    idx = min([length(rawdate) length(data.T1) length(data.T2)]);

    datestr(rawdate(1))
    datestr(cal.time(1))
    CreateFigure;
    [ax, ~] = create_axes(gcf(), 4, 1, 0);

    axes(ax(1)); hold on
    plot(rawdate(1:idx), data.T1(1:idx))
    plot(rawdate(1:idx), data.T2(1:idx))
    legend('T1', 'T2')
    ylabel('T [Volts]')

    axes(ax(2)); hold on
    plot(cal.time, cal.T1)
    plot(cal.time, cal.T2)
    ylabel('T [calibrated]')

    axes(ax(3)); hold on
    plot(data.datenum, data.T1P)
    plot(data.datenum, data.T2P)
    ylabel('Tp [Volts]')

    axes(ax(4)); hold on
    plot(cal.time_tp, cal.T1Pt)
    plot(cal.time_tp, cal.T2Pt)
    ylabel('Tp [calibrated]')
    datetick

    linkaxes(ax, 'x')
    xlabel(ax(4), datestr(cal.time(1), 'mmm-dd'))
    
end