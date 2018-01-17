function [] = pitot_shear(basedir)

    do_pitot = 1;
    do_Tz = 0;

    tfit = 10*60; % in seconds
    min_dz = 0.1; % minimum vertical displacement required for fit

    tempname = [basedir '/proc/temp.mat'];

    if exist(tempname, 'file')
        load(tempname)
    else
        error('temp.mat does not exist. Please run do_temp_proc');
    end

    if ~isfield(T, 'W')
        do_pitot = 0;
    end

    dt = diff(T.time(1:2))*86400;
    if round(dt) > 1
        warning('Time delta in temp.mat is > 1s.')
    end

    N = round(tfit/dt);
    if do_pitot
        head = load('../calib/header_p.mat');
        fits.shear = nan([1, floor(length(T.W)/N)]);
        [spd, ~, ~] = pitot_calibrate(T.W, T.T1, T.P, head.W);
    end
    fits.Tz1 = fits.shear;
    fits.Tz2 = fits.shear;
    fits.time = fits.shear;

    tic;
    for tt=1:N:length(T.W)
        % TODO: support displacements from dp/dt if necessary
        if tt+N > length(T.a_dis_z)
            break;
        end
        zz = -T.a_dis_z(tt:tt+N);
        idx = 1 + (tt-1)/N;

        if std(zz) > min_dz
            if do_pitot
                fits.shear(idx) = robust_regress(zz', spd(tt:tt+N)');
            end
            if do_Tz
                fits.Tz1(idx) = robust_regress(zz', T.T1(tt:tt+N)');
                fits.Tz2(idx) = robust_regress(zz', T.T2(tt:tt+N)');
            end
        else
            if do_pitot
                fits.shear(idx) = nan;
            end
            if do_Tz
                fits.Tz1(idx) = nan;
                fits.Tz2(idx) = nan;
            end
        end
        fits.time(idx) = nanmean(T.time(tt:tt+N));

    end
    toc;

    % save to appropriate files
    if do_pitot
        load ../input/vel_p.mat
        vel_p.shear = interp1(fits.time(1:end-1), fits.shear(1:end-1), vel_p.time)
        save('../input/vel_p.mat', 'vel_p');
    end

end

function [slope] = robust_regress(x,y)

    [coeffs, stats] = robustfit(x-nanmean(x), y-nanmean(y));

    slope = coeffs(2);
end