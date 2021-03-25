function ts_iEEG_writeavg(avg_data,outdir,prefix,varargin)
% function ts_iEEG_writeavg(data_in,outdir,prefix,'option',value,...
%  
% Usage: ts_iEEG_writeavg(avg_data,'./','test','template','tester.avg');
% Writes out a time surfer avg_data structure to an .avg file readable by
% NeuroScan
%
% Required Inputs:
%
% data_in - a valid avg_data structure
% outdir  - directory to write to
% prefix  - prefix of new .avg file(s)
% 
% Optional
%
% template - an template .avg file from same data set
% 
% Based on Information in the following code:
% 
% loadavg - function of eeglab
%
% eeg_write_scan4_avg.m & eeg_load_scan4_avg.m found at
% http://eeg.sourceforge.net/doc_m2html/bioelectromagnetism/
%
% pdf file of head information found at:
% http://dnl.ucsf.edu/users/dweber/dweber_utils/scan_c.html
%
% Created 04/24/2008 by Rajan Patel
% Last Modified on 05/05/2008 by Rajan Patel
%
% See also ts_iEEG_loadavg

%% Check options

if nargin < 3, help(mfilename); end

parms = mmil_args2parms(varargin,...
                        {...
                         'template',[],[],...
                         'logfile',[],[],...
                         'logfid',1,[]...
                       },...
                       false);


errors = ts_checkdata(avg_data);
if ~isempty(errors)
 mmil_error(parms, 'Errors supplied data: %s.', sprintf('\t%s\n', errors{:}));
end

if ~isfield(avg_data,'averages')
    mmil_error(parms,'Input data must be avg_data.');
end

if ~isempty(outdir) && ~exist(outdir,'dir')
 mmil_error(parms,'Directory does not exist: %s.',outdir);
end

for i = 1:length(avg_data.averages)
  outfile{i} = fullfile(outdir,sprintf('%s_%s.avg',prefix,num2str(i)));
  if exist(outfile{i},'file')
    mmil_error(parms,'ERROR: File exists %s.  Change the prefix.',outfile{i});
  end
end

if ~isempty(parms.template) && ~exist(parms.template,'file')
  mmil_error(parms,'Could not find template file: %s.',parms.template);
end

% channel calibration values
SENSITIVITY      = 184.275146; % seems to always be this value
CALIBRATION      = 1; 

% positions in file of certain vars
S_nsweeps_offset 	= 364; % 364 sweep accept (total sweeps 362)
S_pnts_offset      	= 368;
S_nchans_offset 	= 370;
S_variance_offset 	= 375;
S_rate_offset 		= 376;
S_xmin_offset 		= 505;
S_xmax_offset 		= 509;
packed_sizeof_SETUP = 900;
sizeto_DATA         = 900+(75*avg_data.num_sensors);

%% Write File the without a template
if ~isempty(parms.template)
    fid    = fopen(parms.template,'r');
    if fid<0
        mmil_error(parms,'Could not open file for reading: %s.',parms.template);
    end
    fseek(fid,0,'bof');
    header = fread(fid,packed_sizeof_SETUP);
    fseek(fid, S_nchans_offset, 'bof');	  	
    chan = fread(fid, 1, 'ushort');
    fseek(fid,900,'bof');
    for elec = 1:chan
        channel_label_tmp = fread(fid, 10, 'uchar');
        erp1(elec,:) = fread(fid, 47-10, 'uchar')';
        baseline(elec) = fread(fid, 1, 'ushort');
        erp2(elec,:) = fread(fid, 10, 'uchar')';
        sensitivity(elec) = fread(fid, 1, 'float32');
        erp3(elec,:) = fread(fid, 8, 'uchar')';
        calib(elec) = fread(fid, 1, 'float32');
        factor(elec) = calib(elec) * sensitivity(elec) / 204.8;
    end;
    chan_diff = chan-avg_data.num_sensors;
    %% if the number of channels don't match expand the channel info.
    if chan_diff > 0
        for elect = 1:chan_diff
            baseline(end+1)    = baseline(chan);
            sensitivity(end+1) = sensitivity(chan);
            calib(end+1)       = calib(chan);
            factor(end+1)      = factor(chan);
            erp1(end+1,:)        = erp1(chan,:);
            erp2(end+1,:)        = erp2(chan,:);
            erp3(end+1,:)        = erp3(chan,:);
        end       
    end
    channel_info = fread(fid,avg_data.num_sensors*75);
    fclose(fid);
    for i = 1:length(outfile)
        fid = fopen(outfile{i},'w','ieee-le');
        if fid<0
            mmil_error(parms,'Could not open file for writing: %s.',outfile);
        end;
        % CALIBRATION_FACT = avg_data.averages(i).num_trials/(CALIBRATION);
        fseek(fid,0,'bof');
        fwrite(fid,header);
        fseek(fid,S_nsweeps_offset,'bof');
        fwrite(fid,avg_data.averages(i).num_trials,'uint16');
        fseek(fid,S_pnts_offset,'bof');
        fwrite(fid,length(avg_data.averages(i).time),'uint16');
        fseek(fid,S_nchans_offset,'bof');
        fwrite(fid,avg_data.num_sensors,'uint16');
        fseek(fid,S_variance_offset,'bof');
        fwrite(fid,1,'char');
        fseek(fid,S_rate_offset,'bof');
        fwrite(fid,avg_data.sfreq,'uint16');
        fseek(fid,S_xmin_offset,'bof');
        fwrite(fid,avg_data.averages(i).time(1),'float32');
        fseek(fid,S_xmax_offset,'bof');
        fwrite(fid,avg_data.averages(i).time(end),'float32');
        fseek(fid,packed_sizeof_SETUP,'bof');
        fwrite(fid,channel_info);
        for elec = 1:avg_data.num_sensors
          fseek(fid,packed_sizeof_SETUP + (elec-1)*75,'bof');
          tmp_chan           = char(avg_data.sensor_info(elec).label);
          tmp_chan(end+1:10) = 0;
          fwrite(fid,tmp_chan,'uchar');
          fwrite(fid,erp1(elec,:),'uchar');
          fwrite(fid, baseline(elec), 'ushort');
          fwrite(fid,erp2(elec,:),'uchar');
          fwrite(fid, sensitivity(elec), 'float32');
          fwrite(fid,erp3(elec,:),'uchar');
          fwrite(fid,calib(elec), 'float32');
        end
        fseek(fid,packed_sizeof_SETUP + avg_data.num_sensors * 75,'bof');
        for elec = 1:avg_data.num_sensors
            fwrite(fid,zeros(5,1),'char');% unused header for each channel
            fwrite(fid,(avg_data.averages(i).data(elec,:)*avg_data.averages(i).num_trials)/factor(elec),'float32');
        end;
        for elec = 1:avg_data.num_sensors
            fwrite(fid,avg_data.averages(i).stdev(elec,:).^2,'float32');
        end;
        frewind(fid);
        fclose(fid);
    end
else
    %%%%% WITHOUT A TEMPLATE FILE 
    %%%%% AT THIS TIME FSEEK IS USED AS A CHEAT TO MAKE SURE THE KNOWN
    %%%%% INFORMATION IS WRITTEN TO THE RIGHT PLACE
    
    for i = 1:length(outfile)
        fid=fopen(outfile{i},'w','ieee-le');
        if fid<0
            mmil_error(parms,'Could not open file for writing: %s.',outfile);
        end;
        mmil_logstr(parms,'Writing SETUP portion of file.');
        %% WRITE HEADER
        fwrite(fid,zeros(12,0),'uchar');      % rev - length = 12
        fwrite(fid,0,'int32');                % offset to next file
        fwrite(fid,0,'int32');                % offset to prev file
        fwrite(fid,1,'uchar');                % type 1 = avg, 0 = eeg
        fwrite(fid,zeros(20,1),'uchar');      % id - length = 20
        fwrite(fid,zeros(20,1),'uchar');      % oper - length = 20
        fwrite(fid,zeros(20,1),'uchar');      % doctor - length = 20
        fwrite(fid,zeros(20,1),'uchar');      % referral - length = 20
        fwrite(fid,zeros(20,1),'uchar');      % hospital - length = 20
        fwrite(fid,zeros(20,1),'uchar');      % patient - length = 20
        fwrite(fid,0,'int16');                % age
        fwrite(fid,0,'uchar');                % sex
        fwrite(fid,0,'uchar');                % hand
        fwrite(fid,zeros(20,1),'uchar');      % Medications - length = 20
        fwrite(fid,zeros(20,1),'uchar');      % Classification - length = 20
        fwrite(fid,zeros(20,1),'uchar');      % state (Subj Wakefulness)- length = 20
        fwrite(fid,zeros(20,1),'uchar');      % Session label - length = 20
        fwrite(fid,zeros(10,1),'uchar');      % date - length = 10
        fwrite(fid,zeros(12,1),'uchar');      % time - length = 12
        fwrite(fid,0,'float32');              % mean_age (group files only)
        fwrite(fid,0,'float32');              % stdev of age (group files only)
        fwrite(fid,0,'int16');                % n - number in group file
        fwrite(fid,zeros(38,1),'uchar');      % compfile - length = 38 - path and name of comparison file
        fwrite(fid,0,'float32');              % spectwincomp - spectral window compensation factor
        fwrite(fid,0,'float32');              % meanaccuracy - avg response accuracy
        fwrite(fid,0,'float32');              %meanlatency - avg response latency
        fwrite(fid,zeros(46,1),'uchar');      %sortfile - length = 46 (path and name)
        fwrite(fid,0,'int32');                %numevents - number of events in events table
        fwrite(fid,0,'uchar');                %compoper - operation used in comparison
        fwrite(fid,0,'uchar');                %avgmode - set during online averaging
        fwrite(fid,0,'uchar');                %review - set during review of EEG data
        fwrite(fid,0,'uint16');               %nsweeps - number of expected sweeps
        %%% INSERT NUMBER OF TRIALS FROM DATA
        fseek(fid,S_nsweeps_offset,'bof');
        fwrite(fid,avg_data.averages(i).num_trials, 'uint16');%compsweeps - actual number of sweeps
        fwrite(fid,avg_data.averages(i).num_trials, 'uint16');%acceptcnt - accepted number of sweeps
        fwrite(fid,0,'uint16');                               %rejectcnt - number of rejected sweeps
        %%% INSERT NUMBER OF TIME POINTS FROM DATA
        fseek(fid,S_pnts_offset,'bof');
        fwrite(fid,length(avg_data.averages(i).time), 'uint16'); %pnts - number of points per waveform
        %%% INSERT NUMBER OF SENSORS FROM DATA
        fseek(fid,S_nchans_offset,'bof');
        fwrite(fid,avg_data.num_sensors, 'uint16'); %nchannels - number of active channels
        fwrite(fid,0,'uint16');%avgupdate - frequency of avg update
        fwrite(fid,0,'char');%domain - acquisition domain (time = 0, freq = 1)
        fseek(fid,S_variance_offset,'bof');
        fwrite(fid,1, 'char');%variance data included flag
        %% INSERT SAMPLING RATE FROM DATA
        fseek(fid,S_rate_offset,'bof');
        fwrite(fid, avg_data.sfreq, 'uint16'); % sampling rate
        fwrite(fid,0,'double'); %scale factor for calibration
        fwrite(fid,0,'char');%veog corrected flag
        fwrite(fid,0,'char');%heog corrected flag
        fwrite(fid,0,'char');%aux1 corrected flag
        fwrite(fid,0,'char');%aux2 corrected flag
        fwrite(fid,0,'float32');%veog trigger percentage
        fwrite(fid,0,'float32');%heog trigger percentage
        fwrite(fid,0,'float32');%aux1 trigger percentage
        fwrite(fid,0,'float32');%aux2 trigger percentage
        fwrite(fid,0,'int16');%heog chnl number
        fwrite(fid,0,'int16');%veog chnl number
        fwrite(fid,0,'int16');%aux1 chnl number
        fwrite(fid,0,'int16');%aux2 chnl number
        fwrite(fid,0,'char');%veogdir veog trigger direction flag
        fwrite(fid,0,'char');%heogdir heog trigger direction flag
        fwrite(fid,0,'char');%aux1dir aux1 trigger direction flag
        fwrite(fid,0,'char');%aux2dir aux2 trigger direction flag
        fwrite(fid,0,'int16');%veog_n - number points per waveform
        fwrite(fid,0,'int16');%heog_n - number points per waveform
        fwrite(fid,0,'int16');%aux1_n - number points per waveform
        fwrite(fid,0,'int16');%aux2_n - number points per waveform
        fwrite(fid,0,'int16');%veogmaxcnt - number of observations per point
        fwrite(fid,0,'int16');%heogmaxcnt - number of observations per point
        fwrite(fid,0,'int16');%aux1maxcnt - number of observations per point
        fwrite(fid,0,'int16');%aux2maxcnt - number of observations per point
        fwrite(fid,0,'char');%veogmethod - method used to correct
        fwrite(fid,0,'char');%heogmethod - method used to correct
        fwrite(fid,0,'char');%aux1method - method used to correct
        fwrite(fid,0,'char');%aux2method - method used to correct
        fwrite(fid,0,'float32');%ampsensitivity - ext amplifier gain
        fwrite(fid,0,'char');%lowpass - toggle amp low pass filter
        fwrite(fid,0,'char');%highpass - toggle amp high pass filter
        fwrite(fid,0,'char');%notch - toggle amp notch state
        fwrite(fid,0,'char');%autoclipadd - auto add on clipboard
        fwrite(fid,0,'char');%baseline - baseline corrected flag
        fwrite(fid,0,'float32');%offstart - start of baseline correction
        fwrite(fid,0,'float32');%offstop - end of baseline correction
        fwrite(fid,0,'char');%reject - auto reject flag
        fwrite(fid,0,'float32');%rejstart - auto reject start point
        fwrite(fid,0,'float32');%rejstop - auto reject stop point
        fwrite(fid,0,'float32');%rejmin - auto reject min value
        fwrite(fid,0,'float32');%rejmax - auto reject max value
        fwrite(fid,0,'char');%trigtype - trigger type
        fwrite(fid,0,'float32');%trigval - trigger value
        fwrite(fid,0,'char');%trigchnl - trigger channel
        fwrite(fid,0,'int16');%trigmask - wait value for lpt port
        fwrite(fid,0,'float32');%trigisi - interstim interval
        fwrite(fid,0,'float32');%trigmin - min trigger out voltage
        fwrite(fid,0,'float32');%trigmax - max trigger out voltage
        fwrite(fid,0,'char');%trigdir - duration of trigger out pulse
        fwrite(fid,0,'char');%autoscale - autosclae on average
        fwrite(fid,0,'int16');%n2 - number in group 2 (MANOVA)
        fwrite(fid,0,'char');%dir - negative display up or down
        fwrite(fid,0,'float32');%dispmin - y axis
        fwrite(fid,1,'float32');%dispmax - y axis
        %% INSERT FIRST AND LAST TIME POINT FROM DATA
        fseek(fid,S_xmin_offset,'bof');
        fwrite(fid, avg_data.averages(i).time(1), 'float32');   %xmin - epoch start in sec
        fseek(fid,S_xmax_offset,'bof');
        fwrite(fid, avg_data.averages(i).time(end), 'float32'); %xmax - epoch end in sec
        fwrite(fid,0,'float32');%automin - autoscale minimum
        fwrite(fid,0,'float32');%automax - autoscale maximum
        fwrite(fid,0,'float32');%zmin - not used
        fwrite(fid,0,'float32');%zmax - not used
        fwrite(fid,0,'float32');%lowcut - low cut on ext. amp
        fwrite(fid,0,'float32');%highcut - hi cut on ext. amp
        fwrite(fid,0,'char');%common - common mode rejection flag
        fwrite(fid,0,'char');%savemode - eeg, avg or both
        fwrite(fid,0,'char');%manmode - manual rejection of incoming data
        fwrite(fid,zeros(10,1),'char'); %ref - length = 10 - label for reference electrode
        fwrite(fid,0,'char');%rectify - rectification on external channel
        fwrite(fid,0,'float32');%displayxmin - min x axis display
        fwrite(fid,0,'float32');%displayxmax - max x axis display
        fwrite(fid,0,'char');%phase - flag for phase computation
        fwrite(fid,zeros(16,1),'char');%screen - length = 16 - screen overlay path name
        fwrite(fid,0,'int16');%calmode - calibration mode
        fwrite(fid,0,'int16');%calmethod - calibration method
        fwrite(fid,0,'int16');%calupdate - calibration update rate
        fwrite(fid,0,'int16');%calbaseline - baesline correction during calibration
        fwrite(fid,0,'int16');%calsweeps - number of calibration sweeps
        fwrite(fid,0,'float32');%alattenuator - attenuator value for calibr.
        fwrite(fid,0,'float32');%calpulsevolt - voltage for calibr. pulse
        fwrite(fid,0,'float32');%calpulsestart - start time of pulse
        fwrite(fid,0,'float32');%calpulsestop - end time of pulse
        fwrite(fid,0,'float32');%calfreq - sweep frequency
        fwrite(fid,zeros(34,0),'char');%taskfile - length 34 task file name
        fwrite(fid,zeros(34,0),'char');%seqfile - length 34 - sequence file path name
        fwrite(fid,0,'char');%spectmethod - spectral method
        fwrite(fid,0,'char');%spectscaling - scaling employed
        fwrite(fid,0,'char');%spectwindow - window employed
        fwrite(fid,0,'float32');%spectwinlength - length of window
        fwrite(fid,0,'char');%spectorder - order offilter for max entropy method
        fwrite(fid,0,'char');%notchfilter - in or out
        fwrite(fid,0,'int16');%headgain for SYNAMP
        fwrite(fid,0,'int');%additionalfiles - number of add. files
        fwrite(fid,zeros(5,1),'char');%unused - length 5 - free space
        fwrite(fid,0,'int16');%fspstopmethod -FSP stopping method
        fwrite(fid,0,'int16');%fspstopmode - FSP stopping mode
        fwrite(fid,0,'float32');%fspfvalue - FSP F value to term
        fwrite(fid,0,'int16');%fsppoin - FSP single point location
        fwrite(fid,0,'int16');%fspblocksize - FSP block size for averaging
        fwrite(fid,0,'uint16');%fspp1 - FSP start of window
        fwrite(fid,0,'uint16');%fspp2 - FSP end of window
        fwrite(fid,0,'float32');%fspalpha - FSP alpha value
        fwrite(fid,0,'float32');%fspnoise - FSP SNR
        fwrite(fid,0,'int16');%fspv1 - FSP - degrees of freedom
        fwrite(fid,zeros(40,1),'char');%montage - length 40 (file path name)
        fwrite(fid,zeros(40,1),'char');%eventfile - length 40 (file path name)
        fwrite(fid,0,'float32');%fratio - correction factor for spectral array
        fwrite(fid,0,'char');%minor_rev - current minor revision
        fwrite(fid,0,'int16');%eegupdate - how often incomingeeg is refreshed
        fwrite(fid,0,'char');%compressed - data compression flag
        fwrite(fid,0,'float32');%xscale - x pos for scale box (unused)
        fwrite(fid,0,'float32');%yscale - y pos for scale box (unused)
        fwrite(fid,0,'float32');%xsize - waveform size x direction
        fwrite(fid,0,'float32');%ysize - waveform size y direction
        fwrite(fid,0,'char');%acmode - st SYNAP into AC mode
        fwrite(fid,0,'uchar');%commonchnl - channel for common waveform
        fwrite(fid,0,'char');%xtics - scale tool - 'tic' flag in x dir
        fwrite(fid,0,'char');%xrange - scale tool - range(ms,sec,Hz) flag x dir
        fwrite(fid,0,'char');%ytics - scale tool -  'tic' flag in y dir
        fwrite(fid,0,'char');%yrange - scale tool - range (uV,V) flag y dir
        fwrite(fid,0,'float32');%xscalevalue - scale tool - value for x dir
        fwrite(fid,0,'float32');%xscaleinterval - scale tool - interval between tics x dir
        fwrite(fid,0,'float32');%yscalevalue - scale tool - value for y dir
        fwrite(fid,0,'float32');%yscaleinterval - scale tool - interval between tics y dir
        fwrite(fid,0,'float32');%scaletoolx1 - scale tool - upper left hand scrn pos
        fwrite(fid,0,'float32');%scaletooly1 - scale tool - upper left hand scrn pos
        fwrite(fid,0,'float32');%scaletoolx2 - scale tool - lower right hand scrn pos
        fwrite(fid,0,'float32');%scaletooly2 - scale tool - lower right hand scrn pos
        fwrite(fid,0,'int16');%port - port address for ext. triggering
        fwrite(fid,0,'int32');%numsamples - # samples in continuous file
        fwrite(fid,0,'char');%filterflag -
        fwrite(fid,0,'float32');%lowcutoff
        fwrite(fid,0,'int16');%lowpoles - number of poles
        fwrite(fid,0,'float32');%highcutoff
        fwrite(fid,0,'int16');%highpoles
        fwrite(fid,0,'char');%filtertype (bandpass = 0, notch = 1, hipass = 2, lowpass =  3)
        fwrite(fid,0,'char');%filterdomain (frequency = 0, time = 1)
        fwrite(fid,0,'char');%snrflag
        fwrite(fid,0,'char');%coherenceflag
        fwrite(fid,0,'char');%continuoustype - method used to capture events in .cnt
        fwrite(fid,0,'int32');%eventtablepos - position of event table
        fwrite(fid,0,'float32');%continuousseconds - seconds to display per page
        fwrite(fid,0,'int32');%channeloffset - block size of one channel in SYNAMPS
        fwrite(fid,0,'char');%autocorrectflag - auto correct of DC value
        fwrite(fid,0,'uchar');%dcthreshold - auto correct of DC level
        mmil_logstr(parms,'Writing electrode configuration information to file.');
        %% Write electrode info
        for elec = 1:avg_data.num_sensors
            %%% OLD STURCTURE FOR ELECTRODE TABLE
            % char[10] label
            % float32 x_coord
            % float32 y_coord
            % float32 alpha_wt
            % float32 beta_wt
            %% VERSION 3.0 ELECTRODE TABLE
            % write 10 bytes for channel
            %% Insert name from data
            fseek(fid,packed_sizeof_SETUP + (elec-1)*75,'bof');
            tmp_chan           = char(avg_data.sensor_info(elec).label);
            tmp_chan(end+1:10) = 0;
            fwrite(fid,tmp_chan,'uchar');
            fwrite(fid,0,'char');%reference electrode number
            fwrite(fid,0,'char');%skip electrode flag
            fwrite(fid,0,'char');%reject - artifact reject flag
            fwrite(fid,1,'char');%display flag for 'STACK'display
            fwrite(fid,0,'char');%bad electrode flag
            fwrite(fid,avg_data.averages(i).num_trials,'uint16');%n - number of observations
            fwrite(fid,1,'char');%avg_reference status
            fwrite(fid,0,'char');%clipadd - auto add to clipboard
            fwrite(fid,0,'float32');%x_coord for 'TOP' display
            fwrite(fid,0,'float32');%y_coord for 'TOP display
            fwrite(fid,0,'float32');%veog_wt  - VEOG correction weight
            fwrite(fid,0,'float32');%veog_std - VEOG std dev. for weight
            fwrite(fid,0,'float32');%snr statistic
            fwrite(fid,0,'float32');%heog_wt - HEOG correction weight
            fwrite(fid,0,'float32');%heog_std - VEOG std dev for weight
            fwrite(fid,0,'int16');%baseline - baseline correction value in raw ad units
            fwrite(fid,0,'char');%filtered flag
            fwrite(fid,0,'char');%fsp - extra data
            fwrite(fid,0,'float32');%aux1_wt - AUX1 - correction weight
            fwrite(fid,0,'float32');%aux1_std - AUX1 std dev for weight
            fwrite(fid,SENSITIVITY,'float32');%senstivity - of electrode
            fwrite(fid,0,'char');%gain - amplifier
            fwrite(fid,0,'char');%hipass value
            fwrite(fid,0,'char');%lopass value
            fwrite(fid,0,'uchar');%page - display
            fwrite(fid,0,'uchar');%size - electrode window display size
            fwrite(fid,0,'uchar');%impedance test
            fwrite(fid,0,'uchar');%hysicalchnl - physical channel used
            fwrite(fid,0,'char');%rectify - free space
            fwrite(fid,CALIBRATION,'float32');%calibration factor
        end;
        CALIBRATION_FACT = avg_data.averages(i).num_trials/(CALIBRATION);

        % ^---- conversion value to go from microvolts to datapoint
        %
        % convervsion from datapoint to uv is stated as the following:
        % uv = datapoint * calibration / number of sweeps in average
        %
        % however in eeglab when they read in they do the following though:
        % uv = datapoint * (calibration * sensitivity / 204.8) / number of sweeps
        %
        % NOT SURE WHICH ONE! - modify CALIBRATION_FACT if you think it
        % should be changed
        fseek(fid,packed_sizeof_SETUP + avg_data.num_sensors * 75,'bof');
        mmil_logstr(parms,'Writing waveform data to file.');
        for elec = 1:avg_data.num_sensors
            fwrite(fid,zeros(5,1),'char');% unused header for each channel
            fwrite(fid,avg_data.averages(i).data(elec,:)*CALIBRATION_FACT,'float32');
        end;

        for elec = 1:avg_data.num_sensors
            fwrite(fid,avg_data.averages(i).stdev(elec,:).^2,'float32');
        end;

        frewind(fid);
        fclose(fid);

    end
end


