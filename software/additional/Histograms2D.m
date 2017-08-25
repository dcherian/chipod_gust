function Histograms2D(chi, ID, avgfn)

    dt = (chi.time(2)-chi.time(1))*86400;

    CreateFigure;
    ax(1) = subplot(221);
    [counts, xbins, ybins, flag] = hc(chi.dTdz, log10(chi.chi));
    pcolor(xbins, ybins, counts); shading flat;
    xlabel('dT/dz'); ylabel('log_{10} \chi')
    title([ID ' | ' num2str(dt) 's ' avgfn ' estimates'])

    ax(2) = subplot(222);
    [counts, xbins, ybins] = hc(chi.dTdz, log10(chi.Kt));
    pcolor(xbins, ybins, counts); shading flat;
    xlabel('dT/dz'); ylabel('log_{10} K_T')

    ax(3) = subplot(223);
    [counts, xbins, ybins] = hc(chi.dTdz, log10(chi.eps));
    pcolor(xbins, ybins, counts); shading flat;
    xlabel('dT/dz'); ylabel('log_{10} \epsilon')

    ax(4) = subplot(224);
    [counts, xbins, ybins] = hc(chi.dTdz, log10(abs(chi.Jq)));
    pcolor(xbins, ybins, counts); shading flat;
    xlabel('dT/dz'); ylabel('log_{10} |J_q^t|')

    linkaxes(ax, 'x');

    if flag
        title(ax(2), 'Histogram equalized');
    end
    load cmap
    colormap(cmap.chi);
end

function [counts, xbins, ybins, flag] = hc(x, y)

    mask = ~isnan(x) & ~isnan(y);
    [counts, xe, ye] = histcounts2(x(mask), y(mask), ...
                                   'numbins', round(sqrt([numel(x) numel(y)])));

    counts = counts'./max(counts(:));

    % try histogram equalization
    try
        counts = histeq(counts, 200);
        flag = 1;
    catch ME
        flag = 0;
    end

    xbins = (xe(1:end-1) + xe(2:end))/2;
    ybins = (ye(1:end-1) + ye(2:end))/2;
end