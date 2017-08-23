%     this is the driver for Pitot calibration
%
%   created by: 
%        Johannes Becherer
%        Tue Feb 14 11:30:20 PST 2017

clear all;
close all;


%_____________________set processing flags______________________
   do_parallel = 0;     % use paralelle computing 
   do_raw_data = 0;     % do the averaging of the raw-data (1) or skip (0) if done before
   do_v0_self  = 1;     % detremine V0 based on a min of the averaged signal (self contained)
   do_v0_adcp  = 1;     % detremin V0 based on a fit against reference velocity (adcp) data
   do_plot     = 1;     % generate some figures in ../pics/ to compare the different velocity estimates
   do_vel_p    = 0;     % which calibration should be used for vel_p (0 none (default), 1: adcp, 2: self)

%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%_____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name

%_____________________set time limits______________________
 % get time limits from whoAmI;
   [TL] =   whoAmI_timeLimits(basedir);
   time_range      = TL.pitot;
 % set manually
   %time_range(1)  = datenum(2000, 1, 1, 0, 0, 0);
   %time_range(2)  = datenum(2030, 1, 1, 0, 0, 0);

 % This is the time range where the pitot sensor is returning
   % good data
   time_range(1)  = datenum(2014, 12, 04, 14, 0, 0);
   time_range(2)  = datenum(2015, 11, 15, 15, 0, 0);

   % calibrate in time range different from valid data time range?
   % if so set limits here just as for time_range.
   % by default, both time ranges are equal.
   cal_time_range = time_range;

   % which temperature sensor to use T1 (1) or if T1 is broken T2 (2) ;  
   % for gusTs (0)
   
   if isfield(TL, 'T1')     % chipod
      % check which sensor longer survives
      if TL.T1(2) < TL.pitot(2) & TL.T2(2)>TL.T1(2)
         use_T = 2;
         % adjust pitot calibration time if temp dies early
         if TL.T2(2) < cal_time_range(2) 
            cal_time_range(2) =  TL.T2(2);
         end
      else
         use_T = 1;
         % adjust pitot calibration time if temp dies early
         if TL.T1(2) < cal_time_range(2) 
            cal_time_range(2) =  TL.T1(2);
         end
      end
   
   else  % gust
     use_T = 0; 
      % adjust pitot calibration time if temp dies early
      if TL.T(2) < cal_time_range(2) 
         cal_time_range(2) =  TL.T(2);
      end
   end

   % manually
   % use_T = 1;


   % shall the pressure calibration be switched off

      if TL.P(2) == TL.master(1)
         use_press = 0;
      else
         use_press = 1;
      end
      % set flag manually 
      %use_press   =  0;
   
      % adjust pitot calibration time if pressure dies early
      if TL.P(2) < cal_time_range(2) & use_press
         cal_time_range(2) =  TL.P(2);
      end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% DO NOT CHANGE BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%_____________________do raw data processing______________________
   if do_raw_data
      generate_praw( basedir, do_parallel, time_range);
   end

   % % pressure
   % P.P    =  Praw.P*head.coef.P(2) + head.coef.P(1);

   % % compass
   %    if isfield(head.coef, 'CMP')
   %       P.cmp  = Praw.cmp + head.coef.CMP(1);
   %    else
   %       P.cmp  = Praw.cmp;
   %       disp(['CMP' ' does not exit in header']);
   %    end

   % %---------------------pre calibration for Pitot----------------------
   %    %% find all idexes in the desired time interval;
   %    iiPcal = find( P.time>=cal_time_range(1) & P.time<=cal_time_range(2) );
   %    iiP = find( P.time>=time_range(1) & P.time<=time_range(2) );

   % % set the average temperature as reference value for the Pitot calibration
   % W.T0   =  nanmean(P.T(iiPcal));
   % W.P0   =  nanmean(P.P(iiPcal));

   % % calibrate the Pitot voltage for temperature (pressure ? Tilt ?)
   % P.W   =   Praw.W - (P.T-W.T0)*W.T(2);

   % % take out long-term drifts.
   % plot(P.time, P.W);
   % % [t, v] = ginput()
   % tpts = [7.359373271889401
   % 7.359437788018433
   % 7.359907834101382
   % 7.360156682027649
   % 7.360359447004609
   % 7.360433179723503
   % 7.362165898617512
   % 7.362267281105991
   % 7.362470046082949
   % 7.362829493087558] * 1e5;

   % for tt=1:length(tpts)
   %     it = find_approx(P.time, tpts(tt));
   %     Ppts(tt) = P.W(it);
   % end
   % hold on; plot(tpts, Ppts, 'k--')

%_____________________determine V0______________________
   
   if do_v0_self | do_v0_adcp | do_plot
      determine_v0( basedir, do_v0_self, do_v0_adcp, do_plot, do_vel_p, time_range, use_T, use_press )
   end

 
