%     this is the driver for the prep-processing for  chi-processing
%
%   created by: 
%        Johannes Becherer
%        Tue Sep 20 10:51:19 PDT 2016

clear all;
close all;


%_____________________set processing flags______________________
   do_parallel = 0;     % use paralelle computing 
   do_temp     = 0;     % generate temp.mat 
   do_vel_p    = 0;     % generate vel_p.mat
   do_vel_m    = 0;     % generate vel_m.mat
   do_dTdz_m   = 1;     % generate dTdz_m.mat
   do_dTdz_i   = 0;     % generate dTdz_i.mat 
   use_pmel    = 0;     % use TAO/TRITON/PIRATA/RAMA mooring data?
   use_rama    = 1;     % use prelim processed RAMA data
   RamaPrelimSalCutoff = 1/(1*60*60); % filter cutoff (Hz) for
                                      % filtering prelim RAMA
                                      % salinity data (set NaN to disable)
   use_TS_relation = 0; % fit TS relation to estimate N2 from
                        % mooring data?
   use_mooring_sal = 1; % use mooring salinity along with dTdz_i
                        % to estimate N^2 in dTdz_i.
                        % otherwise code assumes fixed salinity=35.
   modify_header_info = 0; % fix header info in this file?

%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name

%_____________________get list of all raw data______________________
   [fids, fdate] = chi_find_rawfiles(basedir);

%_____________________for automated PMEL mooring processing____________
   % chipod location (positive North, East & Down)
   ChipodLon = 90; ChipodLat = 12; ChipodDepth = 15;

   if use_pmel
       pmeldir = '~/TaoTritonPirataRama/'; % directory with pmel mooring files
                                            % (can obtain an updated copy from ganges)
       % which high-freq data file should I use?
       % 2m/10m/30m/hr
       velfreq = '30m';
       Tfreq = '10m';
       Sfreq = 'dy';

       % find start and end of depoyment from raw files
       rawdir = [basedir filesep 'raw' filesep];
       data = raw_load_chipod([rawdir fids{1}]);
       deployStart = data.datenum(1);
       data = raw_load_chipod([rawdir fids{end}]);
       deployEnd = data.datenum(end);
   end

   % use RAMA preliminary data
   if use_rama
       ramaname = '~/rama/RamaPrelimProcessed/RAMA13.mat';
   end

%%%%%%%%%%%%%%% header data %%%%%%%%%%%%%%%%%%%%%%

if modify_header_info
    fid = [basedir filesep 'calib' filesep 'header.mat'] ;

    % read raw-data
    [~, head] = raw_load_chipod([rawdir fids{1}]);

    % head.irep values from raw data appear to be crap.
    % delete them so that I'm sure they're not being used.
    % johannes says the RAMA deployment is on standardized smapling
    % frequencies, so the code automagically makes the right
    % choice.
    head.irep = [];

    % fix compass offset + declination
    % this is *only* used in chi_calibrate_chipod.m
    % Compass is calibrated to true N in Corvallis
    % Corvallis declination is 15° 44' °E in 2013.
    % (not sure that this is relevant).
    % https://www.ngdc.noaa.gov/geomag-web/#declination
    % RAMA13 (SJW) declination at 12N, 90E is 1.12°W
    head.coef.CMP(1) = -101.15  + 1.12;

    % save header in proper destination
    save(fid, 'head');

    %________________ pitot calibrations __________________
    fid = [basedir filesep 'calib' filesep 'header_p.mat'] ;
    W.T  = [0 -0.003154669 0 0 0];
    W.Ps = [0 0 0 0 0]; % pressure sensor is bad; accounted for in offset
    W.Tilt = [0 0.000088684 0 0 0];
    W.Pd = [0 0.0003995 0 0 0]; %if slope>1 else W.Pd = [0 slope 0 0 0];
    assert(W.Pd(2) < 1, 'WPd(2) > 1 !');

    % offsets
    W.V0 = 0;
    W.P0 = 0;
    W.T0 = 0;

    save(fid, 'W');

    % % below is from sally's cali_chipod_rama13.m
    % % saved here as backup
    % % in this deployment, a lot of the header is incorrect.
    % the data is recorded at 100 Hz for time and T1P and T2P
    % the data is recorded at 50 Hz for T1, T2, AX, AY, AZ, W, WP, P
    % the data is recorded at 5 Hz for CMP
    % head.irep.T1P   = 20;
    % head.irep.T1    = 10;
    % head.irep.AX    = 10;
    % head.irep.AY    = 10;
    % head.irep.AZ    = 10;
    % head.irep.W1    = 10;
    % head.irep.W     = 10;
    % head.irep.WP    = 10;
    % head.irep.T2P   = 20;
    % head.irep.T2    = 10;
    % head.irep.P     = 10;
    % head.irep.R2    = 10;
    % head.irep.R3    = 10;
    % head.irep.R4    = 1;
    % head.irep.R1    = 1;
    % head.irep.VA    = 1;
    % head.irep.MK0   = 1;
    % head.irep.MK1   = 1;
    % head.irep.QUE   = 1;
    % head.irep.VD    = 1;
    % head.irep.CMP   = 1;
    % head.irep.MK5   = 1;
    % head.irep.MK6   = 1;
    % head.oversample and samplerate are wrong too, but all the variables that
    % we need are correct.
    % head.oversample = [1 2 2 2 2 2 2 2 1 2 2 2 2 20 20 20 20 20 20 20 20 20 20];
    % head.samplerate = [100 50 50 50 50 50 50 50 100 50 50 50 50 5 5 5 5 5 5 5 5 5 5];
end

%%%%%%%%%%%%%%%%%%% temp processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
if do_temp
   %_____________processing loop through all raw files__________________

      % init parallel pool
      if(do_parallel)
         parpool;
         % parallel for-loop
         parfor f=1:length(fids)
            try % take care if script crashes that the parpoo is shut down
               disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
               chi_T_proc(basedir, fids{f});
            catch
               disp(['!!!!!! ' fids{f} ' crashed while processing T structure !!!!!!' ]);
            end
         end
         % close parpool
         delete(gcp);
      else
         for f=1:length(fids)
            disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
            chi_T_proc(basedir, fids{f});
         end
      end

   %____________________merge individual files______________________
      
      % average 20 sec
      chi_merge_and_avg(basedir, 'temp', 20);
end


%%%%%%%%%%%%%%%%%%% generating Pitot velocity input file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
if do_vel_p

   fidf = '../proc/P_fit.mat';
   fids = '../proc/P_self.mat';
   
   if exist(fidf, 'file');
      load(fidf);
      vel_p.text = 'vel_p.mat is generated based on the ADCP fitted Pitot signal';
      disp(vel_p.text);
   elseif exist(fids, 'file');
      load(fids);
      vel_p.text = 'vel_p.mat is generated in the self contained way';
      disp(vel_p.text);
   else
      disp([fid ' does not exist. Run calibrate_pitot first !']);
   end

   vel_p.time  = P.time;
   vel_p.spd   = P.spd;
   vel_p.U     = P.U;
   vel_p.u     = real(P.U);
   vel_p.v     = imag(P.U);

   save('../input/vel_p.mat', 'vel_p');
   
end


%%%%%%%%%%%%%%%%%%% mooring velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
if do_vel_m
    sdir  = [basedir filesep 'input' filesep];

    if use_pmel
        moor = ExtractUVFromTaoTritonPirataRama(ChipodLon, ChipodLat, ...
                                                ChipodDepth, deployStart, ...
                                                deployEnd, pmeldir, ...
                                                'RAMA', velfreq);
    end

    %_______ EXAMPLE________________
    % load('../../../mooring_data/mooring_Pirata14_524.mat') ;

    chi_generate_vel_adcp(moor.time, moor.depth, moor.u, moor.v, moor.depth, sdir);
end


%%%%%%%%%%%%%%%%%%% mooring dTdz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
if do_dTdz_m
      sdir  = [basedir filesep 'input' filesep];

      if use_pmel
          [T1, T2] = ExtractTSFromTaoTritonPirataRama(ChipodLon, ChipodLat, ...
                                                      ChipodDepth, deployStart, ...
                                                      deployEnd, pmeldir, 'RAMA', ...
                                                      Tfreq, Sfreq);
          save([basedir filesep 'proc' filesep 'T_m.mat'], ...
                'T1', 'T2')
      end

      if use_rama
          % Using RAMA prelminary data (10 min T,S)
          [T1, T2] = ExtractTSFromRamaPrelim(ramaname, ChipodDepth);

          if ~isnan(RamaPrelimSalCutoff)
              disp('Low pass filtering RAMA salinity')
              Sfilt = gappy_filt(1./diff(T1.time(1:2)*86400), ...
                                 {['l' num2str(RamaPrelimSalCutoff)]}, ...
                                 4, T1.S);
              T1.S = Sfilt;
              % figure; plot(T1.time, T1.S); hold on; plot(T1.time, Sfilt);

              Sfilt = gappy_filt(1./diff(T2.time(1:2)*86400), ...
                                 {['l' num2str(RamaPrelimSalCutoff)]}, ...
                                 4, T2.S);
              T2.S = Sfilt;
          end
      end

      %_______ EXAMPLE________________
      %  load('../../G002/proc/temp.mat') ; % surounding instruments
      %     T1.time = T.time; 
      %     T1.z    = nanmedian(T.depth); 
      %     T1.T    = T.T; 
      %     T1.S    = ones(size((T.T)))*35; 
      %  load('../../G011/proc/temp.mat') ; % surounding instruments
      %     T2.time = T.time; 
      %     T2.z    = nanmedian(T.depth); 
      %     T2.T    = T.T; 
      %     T2.S    = ones(size((T.T)))*35; 

      chi_generate_dTdz_m(T1.time, T1.z, T1.T, T1.S, ...
                          T2.time, T2.z, T2.T, T2.S, sdir, use_TS_relation);

      %__________________recalculate N^2 using processed mooring salinity____________________

      if use_mooring_sal
          if ~exist('../input/dTdz_i.mat', 'file')
              error(['Create dTdz_i.mat first. Run pre_driver with ' ...
                     'do_dTdz_i=1']);
          end

          load ../input/dTdz_i.mat
          load ../input/dTdz_m.mat

          % interpolate to Tz_i.time
          dSdz = interp1(T1.time, (T1.S-T2.S)/abs(T1.z-T2.z), Tz_i.time);
          Smean = interp1(T1.time, (T1.S + T2.S)/2, Tz_i.time);

          Tnames = {'T1', 'T2', 'T12'};
          Tznames = {'Tz1', 'Tz2', 'Tz12'};
          Nnames = {'N2_1', 'N2_2', 'N2_12'};
          for ii=1:length(Tnames)
              Tii = Tz_i.(Tnames{ii});
              dTdz = Tz_i.(Tznames{ii});

              alpha = sw_alpha(Smean, Tii, ChipodDepth);
              beta = sw_beta(Smean, Tii, ChipodDepth);

              Tz_i.(Nnames{ii}) = -9.81 * (-alpha.*dTdz + beta.*dSdz);
          end

          Tz_i.S = Smean;
          save('../input/dTdz_i.mat', 'Tz_i');
      end
end

%%%%%%%%%%%%%%%%%%% internal dTdz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
if do_dTdz_i
   %_____________processing loop through all raw files__________________

      dt   = 60; % sec bits of data for analysis
      do_P = 0; % use pressure instead of acceleration to get z 

      disp('calculating the intrenal dTdz');
      % init parallel pool
      if(do_parallel)
         parpool;
         % parallel for-loop
         parfor f=1:length(fids)
            try % take care if script crashes that the parpoo is shut down
               disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
               chi_generate_dTdz_i(basedir, fids{f}, dt, do_P);
            catch
               disp(['!!!!!! ' fids{f} ' crashed while processing  internal dTdz structure !!!!!!' ]);
            end
         end
         % close parpool
         delete(gcp);
      else
         for f=1:length(fids)
            disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
            chi_generate_dTdz_i(basedir, fids{f}, dt, do_P);
         end
      end

   %____________________merge individual files______________________
      chi_merge_and_avg(basedir, 'dTdz', 600);

   %_____________________cp result to the input directory______________________
   ! cp ../proc/dTdz.mat ../input/dTdz_i.mat
end
