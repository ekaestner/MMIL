%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Align Intan micromacro data with clinical sEEG data
%
% Charlie Dickey
% cdickey@ucsd.edu
% 04/09/2019
% Halgren Lab
% 
% This script was written to help analyze data from SD018. This patient had
% both sEEG and micromacro data recorded, but the triggers were only 
% recorded in the sEEG clinical system.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% load sEEG data from micromacro electrode and trigger channel

% SEEG sampling rate
fs_clin = 1024;

% sEEG file containing HH task data 
fname = '/space/seh8/3/halgdev/projects/bqrosen/MULTI_initial/seeg_UCSD/SD018/SD018_Day07.edf';

% optional: load only the header to inspect metadata
header = edfread(fname); 

% load referential sEEG data recorded from the macro contacts of the
% micromacroelectrode (Macro1-7) as well as the trigger data (TRIG)
targetSignals = {'Macro1', 'Macro2', 'Macro3', 'Macro4', 'Macro5', 'Macro6', 'Macro7', 'TRIG' };
tic
[header, data_clin] = edfread(fname, 'targetSignals', targetSignals);
toc

% isolate trigger and sEEG data
data_trig = data_clin(8,:);
data_clin = data_clin(1:7,:);

% apply 60 Hz notch filter with harmonics up to Nyquist
for nf = [60,120,180,240,300,360,420,480]
    Wo = nf/(fs_clin/2);
    BW = Wo/35;
    [b,a] = iirnotch(Wo, BW);
    data_clin = filtfilt(b,a,data_clin')';
end

% compute whole probe bipolar
data_clin_wp = data_clin(1,:)-data_clin(7,:);


%% load bipolar Intan micromacro data

% Intan file directory
intan_dir = '/space/seh9/1/halgdev/projects/micromacro/3_24_19';

% Intan files (each 1 min long) containing HH task data
intan_files = {'SD018_Sun1_uM__190324_132409.rhd', 'SD018_Sun1_uM__190324_132509.rhd', 'SD018_Sun1_uM__190324_132609.rhd',...
    'SD018_Sun1_uM__190324_132709.rhd', 'SD018_Sun1_uM__190324_132809.rhd', 'SD018_Sun1_uM__190324_132909.rhd',...
    'SD018_Sun1_uM__190324_133009.rhd', 'SD018_Sun1_uM__190324_133109.rhd', 'SD018_Sun1_uM__190324_133209.rhd',...
    'SD018_Sun1_uM__190324_133309.rhd'};

% load and concatenate Intan files from list
data_intan = [];
fs_intan = [];
for file = 1:numel(intan_files)
    disp(string(file) + '/' + string(length(intan_files)) + ' complete')
    data_intan_tmp = read_Intan_RHD2000_file_nogui(sprintf('%s/%s', intan_dir, intan_files{file})); % day 4
    data_intan = [data_intan data_intan_tmp.amplifier_data];
    fs_intan = [fs_intan data_intan_tmp.frequency_parameters.amplifier_sample_rate];
end

clear data_intan_tmp

% check whether sample rate for all Intan files is the same
if range(fs_intan) == 0
    fs_intan = fs_intan(1);
else
    disp('ERROR: Sample rates differ across Intan files.')
end

% 60 Hz notch filter with harmonics up to 480
for nf = [60,120,180,240,300,360,420,480]
    Wo = nf/(fs_intan/2);
    BW = Wo/35;
    [b,a] = iirnotch(Wo, BW);
    data_intan = filtfilt(b,a,data_intan')';
end


%% align data using whole probe bipolar

% channel 16 = whole probe bipolar (Macro1 - Macro7)
intan_ch = 16;

% downsample Intan data to clinical sampling rate (includes auto anti-aliasing)
data_intan_ds = resample(data_intan', fs_clin, fs_intan)';

% compute cross-correlation and find alignment time lag based on max value
xc = xcorr(data_clin_wp, data_intan_ds(intan_ch,:));
[~,lag] = max(xc);
offset = lag-size(data_clin_wp,2)+1; % use this value for the alignment


%% visualize alignment

% window to visualize
start = 3*fs_clin;
stop = 5*fs_clin;
times = (1:stop-start+1)./fs_clin;

% plot alignment
hold on
plot(times, data_intan_ds(intan_ch, start:stop), 'r')
plot(times, data_clin_wp(start+offset:stop+offset), 'b')
xlim([0 times(end)])
title('alignment of Intan and clinical recordings')
legend('intan whole probe bipolar', 'clinical whole probe bipolar', 'Location', 'southwest')
xlabel('time (s)')
ylabel('amplitude (\muV)')
hold off
