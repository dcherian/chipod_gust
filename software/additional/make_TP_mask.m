function [mask] = make_TP_mask(T, time_range)

    time_range(1) = datenum(2014, 10, 8);
    time_range(2) = datenum(2014, 10, 14);

    tt = [find_approx(T.time, time_range(1)): ...
          find_approx(T.time, time_range(2))];

    plot(T.time(tt), T.T1Pt(tt));

    datetick;
end
