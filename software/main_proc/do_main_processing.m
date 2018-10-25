function [] = do_main_processing( basedir, pflag, time_lim, do_just_merge )
%% [] = do_main_processing( basedir, [pflag], [time_lim], [do_just_merge] ) 
%  
%     This function drives the main processing for a given GusT or Chipod.
%
%     INPUT
%        basedir        :  base directory of the given instrument
%        pflag          :  processing flag object generated by chi_processing_flags.m
%                          if pflag is not given it will be generated automatically based on 
%                          inputs in the base directory
%        time_lim       :  time limits (notes this effects only the merged product not the processing time)
%        do_just_merge  :  shall the processing be skiped to do just the merging? (default 0) 
%
%   created by: 
%        Johannes Becherer
%        Thu Dec 28 12:33:57 PST 2017

ticstart = tic;

if nargin < 4
   do_just_merge = 0;
end


if nargin < 3 % if no time limits are set we use pratical no limits
   time_lim = [datenum(1900,1,1) datenum(2100,1,1)];
end

%____________________set automatic pflag______________________
   if nargin < 2

      pflag = chi_processing_flags;     % get list of processing flags

      %---------------set processing flags automatically----------
      pflag = pflag.auto_set(basedir);
      pflag.master.parallel = 0;
      %---------------------get flag status----------------------
      pflag.status();

   end

%_____________________do main processing______________________
   %_____________________get all raw files______________________

   if ~do_just_merge
      [fids, fdate] = chi_find_rawfiles(basedir, time_lim);

      if(pflag.master.parallel)

         % open the parallel pool
         if pflag.master.parallel == 1
            parpool;
         else
            parpool(pflag.master.parallel);
         end

         % parallel for-loop
         parfor f=1:length(fids)
             disp(['processing day ' num2str(f) ' of ' num2str(length(fids))]);
             try % take care if script crashes that the parpoo is shut 
                filedate = datenum(fdate{f}(1:8), 'yymmddhh');
                if (filedate) >= floor(time_lim(1)) & (filedate) < ceil(time_lim(2))
                  chi_main_proc(basedir, fids{f}, pflag);
                end
             catch ME
                 disp(['!!!!!! ' fids{f} ' (f = ' num2str(f) ') crashed while ' ...
                       'processing  !!!!!! \n\n' ME.message]);
             end
         end
         % close parpool
         delete(gcp);
      else
          for f=120:length(fids)
             disp(['processing day ' num2str(f) ' of ' num2str(length(fids))]);
             try
                filedate = datenum(fdate{f}(1:8), 'yymmddhh');
                if (filedate) >= floor(time_lim(1)) & (filedate) < ceil(time_lim(2))
                      chi_main_proc(basedir, fids{f}, pflag);
                 end
             catch ME
                 disp(['!!!!!! ' fids{f} ' (f = ' num2str(f) ') crashed while ' ...
                       'processing  !!!!!! \n\n' ME.message]);
             end
         end
      end
   end

   %_____________________merge all days______________________
   disp('merge all days and estimates.')
   ticmerge = tic;
      %_loop through all processing flags for chi processing
      % keep averaging window 0 here.
      % Only merge, average later in combine_turbulene.m
      for i = 1:length(pflag.id)
            [id, ~, ~, ~] = pflag.get_id(i);
            if pflag.proc.(id) % check if flag is active
               ddir = ['chi' filesep 'chi_' id];
               chi_merge_and_avg(basedir, ddir, 0, time_lim );

               if ~contains(id, '_ic')
                   try
                       ddir = ['chi' filesep 'chi_' id filesep 'stats'];
                       mkdir([basedir filesep 'proc' filesep 'chi' filesep 'stats' filesep]);
                       statsname = ['chi' filesep 'stats/chi_' id '_stats'];
                       chi_merge_and_avg(basedir, ddir, 0, [], statsname);
                   catch ME
                       disp(ME)
                       disp('Error! have the fitting stats been saved? Skipping...')
                   end
               end
            end
      end
   disp('Finished merging all estimates.')
   toc(ticmerge);

   %_____________merge eps data______________________
   if pflag.master.epsp
      % keep averaging window 0 here.
      % Only merge, average later in combine_turbulene.m
      chi_merge_and_avg(basedir, 'eps', 0, time_lim);
   end

   disp('Finished running main processing.')
   toc(ticstart)
