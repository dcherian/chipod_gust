addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%____________________set directories______________________
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name
   rawdir       = [basedir filesep 'raw' filesep]; % raw files location

% will create header.mat file if necessary
% if header.mat exists, it will read it
head = chi_get_calibration_coefs(basedir);

% will create header.mat file if necessary
% if header_p.mat exists, it will read it
try
    W  = chi_get_calibration_coefs_pitot(basedir);
catch ME
    warning('Could not load pitot calibration')
end

format long;
disp(' ')
disp('----- temp T1; T2')
[head.coef.T1'; head.coef.T2']

disp(' ')
disp('----- differentiator T1P; T2P')
[head.coef.T1P'; head.coef.T2P']

disp(' ')
disp('----- Accelerometer: AX; AY: AZ')
[head.coef.AX'; head.coef.AY'; head.coef.AZ']

disp(' ')
disp('----- Compass')
disp('First element might have been corrected for declination.')
head.coef.CMP'


disp(' ')
disp('----- Pressure')
head.coef.P'

if exist('W', 'var')
    disp(' ')
    disp('----- Pitot')
    W
end

format short;