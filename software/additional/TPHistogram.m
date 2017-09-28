function [] = TPHistogram(TP, Trange, mask, masklabel)

    if ~exist('Trange', 'var'), Trange = 1:length(TP.time); end
    if ~exist('masklabel', 'var'), masklabel = 'masked'; end

    CreateFigure;

    subplot(121);
    histogram(log10(TP.T1Pvar(Trange)), 'FaceColor', [1 1 1]*0.8, ...
              'displayname', 'all')
    hold on;
    if exist('mask', 'var')
        tp = TP.T1Pvar(Trange);
        histogram(log10(tp(mask)), 'displayname', masklabel)
    end
    xlabel('log_{10} var(T1P)')
    ylabel('count')
    legend('-dynamiclegend');

    subplot(122)
    histogram(log10(TP.T2Pvar(Trange)), 'FaceColor', [1 1 1]*0.8, ...
              'displayname', 'raw')
    hold on;
    if exist('mask', 'var')
        tp = TP.T2Pvar(Trange);
        histogram(log10(tp(mask)), 'displayname', 'masklabel')
    end
    xlabel('log_{10} var(T2P)')
    ylabel('count')