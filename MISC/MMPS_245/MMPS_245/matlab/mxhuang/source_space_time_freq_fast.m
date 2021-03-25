% source_space_time_freq_fast
% This program perform source-space wavelet-based
% time and frequency analysis for single trial evoked responses
% 
% before running this program, you should run "dSPM_bem_vectorview"
% on average response to obtain the inverse operator
%
% (c) Mingxiong Huang, Ph.D. 03/23/2006

clear;                                  	% clear working space
INPUT = 'source_space_time_freq_default_input';	% file to get default inputs
TEMP = 'source_space_time_freq_temp_output';	% file to write this run's inputs
set(0,'DefaulttextInterpreter','none');

wd = pwd;			% get current working directory

if(exist([wd '/' INPUT]) == 2), % a defaults file exists from a prior run
  f_default = fopen(INPUT,'r');
else
  f_default = [];		% none found in current directory
end

if(exist([wd '/' TEMP]) == 2),
  delete temp_output		% we'll write new inputs to here
end
f_temp = fopen(TEMP,'w');


mn_inverse_file=default_s('Name of dSPM result file (.mat)',f_default,f_temp);
eval(['load ',mn_inverse_file,' id_grad mn_inverse id_lh_dip id_rh_dip n_grid']);

raw_filename=default_s('Name of single trial MEG file (.fif)',f_default,f_temp);

result_file=default_s('Name of result file (no ext.)',f_default,f_temp);

stc_lh_file_Z=['stc/',result_file,'_Z-lh.stc'];
stc_rh_file_Z=['stc/',result_file,'_Z-rh.stc'];
stc_lh_file_F=['stc/',result_file,'_F-lh.stc'];
stc_rh_file_F=['stc/',result_file,'_F-rh.stc'];

t_interval=default_v('Pre-, post-stim interval (e.g. [-1000 2000])',f_default,f_temp);
tpre=t_interval(1);
tpost=t_interval(2);
t_base=default_v('Baseline correction interval(e.g. [-800 0])',f_default,f_temp);
t_save=default_v('Interval to save (e.g. [-800 1500])',f_default,f_temp);
f0=default_v('Frequency of interest (e.g. 4)',f_default,f_temp);
width=default_v('Wavelet width (e.g. 1)',f_default,f_temp); % sigma_f=f0/width
EOG_th=default_v('EOG threshold (e.g. 150)',f_default,f_temp);
MEG_grad_th=default_v('MEG Gradiometer threshold (e.g. 3000)',f_default,f_temp);
cond_avg=default_v('The number of condition to average (e.g. [1 2 5 8] for all new, [4 6 16 32] for all old)',f_default,f_temp);

%**************** STEP ***************************************************
% Step:  We're done with user input.  Save responses for next time

if(~isempty(f_default)),	% there was an input file
  fclose(f_default);
  eval(['delete ' INPUT]);	% remove it
end
fclose(f_temp);
computer_type=computer;
if strcmp(computer_type,'PCWIN')==1
   eval(['!copy ' TEMP ' ' INPUT ]);
else   
	eval(['! mv ' TEMP ' ' INPUT ]); % move to input for next time
end

% detect the triggers
fprintf('Detecting the trigers...');
ev2_list=ez2read(raw_filename);
ev2_list(:,2)=ev2_list(:,2)-1; % being consistent with Mao's setting
fprintf('DONE\n')

% load the raw data
fprintf('Reading the raw single trial data...');
B=read_raw_fif(raw_filename);
fprintf('DONE\n')

%cond_avg=[1 2 5 8]; % the number of conditions to average, e.g. all new
%cond_avg=[4 6 16 32]; % all old conditions
id_all_avg=[];
for i=1:length(cond_avg)
    id_all_avg=[id_all_avg;find(ev2_list(:,1)==cond_avg(i))];
end    

sf=rawdata('sf',raw_filename);

np_pre=round(tpre/1000*sf);
np_post=round(tpost/1000*sf);

% reshape B
B=reshape(B,size(B,1),size(B,2)*size(B,3));

% get channel names
[NA,KI,NU,CT,T]=channames(raw_filename);

id_EOG=find(KI==202); % find the EOG channel

% now loop through different single trials
fprintf('Now loop through the trigers for time-freq analysis...\n');
trial_good=[];
time_pre_post=((np_pre):np_post)*1000/sf;
index_base=find(time_pre_post>t_base(1) & time_pre_post < t_base(2));
index_save=find(time_pre_post>=t_save(1) & time_pre_post <=t_save(2));
TFR_avg=single(zeros(size(mn_inverse,1),length(time_pre_post)));
TFR_sensor=zeros(length(id_grad),length(time_pre_post));
B_avg=zeros(length(id_grad),length(time_pre_post));
J_trial=TFR_avg;
h_waitbar = waitbar(0,'Single trial time-freq analysis, please wait...');
for i=1:length(id_all_avg)
    fprintf('Trial %d of total %d trials\n',i,length(id_all_avg));
    id_trial_pre_post=(ev2_list(id_all_avg(i),2)+np_pre):(ev2_list(id_all_avg(i),2)+np_post);
    B_trial=B(:,id_trial_pre_post);
    B_baseline=mean(B_trial(:,index_base)')';
    B_meg=B_trial(id_grad,:)-B_baseline(id_grad)*ones(1,length(id_trial_pre_post)); % baseline correction
    if (max(abs(B_trial(id_EOG,:)))*1e6 < EOG_th) & (max(abs(B_meg(:)))*1e13 < MEG_grad_th) % reject artifacts
        trial_good=[trial_good;i];
        TF_B_meg=zeros(size(B_meg));
        for j=1:size(B_meg,1)
            TF_B_meg(j,:)=traces2TF_complex(B_meg(j,:),f0,sf,width); % calculate the wavelet convolution
        end
        TFR_sensor=TFR_sensor+abs(TF_B_meg).^2;
        B_avg=B_avg+B_meg;
        %TF_B_baseline=mean(TF_B_meg(:,index_base)')';
        %TF_B_meg=TF_B_meg-TF_B_baseline*ones(1,size(TF_B_meg,2)); % baseline correction
        J_trial=mn_inverse*single(TF_B_meg); % the inverse to source space    
        TFR_avg=TFR_avg+abs(J_trial).^2;
    end
    waitbar(i/length(id_all_avg),h_waitbar);
end
TFR_avg=TFR_avg/length(trial_good);
TFR_sensor=TFR_sensor/length(trial_good);
B_avg=B_avg/length(trial_good);
close(h_waitbar)

clear J_trial TFR_trial B mn_inverse

% write the stc files
TFR_avg_all_in_one=TFR_avg(1:3:size(TFR_avg,1),:)+TFR_avg(2:3:size(TFR_avg,1),:)+TFR_avg(2:3:size(TFR_avg,1),:);
TFR_avg_all_in_one_base=mean(TFR_avg_all_in_one(:,index_base)')';
TFR_avg_all_in_one_base_std=std(TFR_avg_all_in_one(:,index_base)')';
TFR_avg_all_in_one_Z=TFR_avg_all_in_one-TFR_avg_all_in_one_base*ones(1,size(TFR_avg_all_in_one,2));
TFR_avg_all_in_one_Z=TFR_avg_all_in_one_Z./(TFR_avg_all_in_one_base_std*ones(1,size(TFR_avg_all_in_one,2)));
TFR_avg_all_in_one_F=TFR_avg_all_in_one./(TFR_avg_all_in_one_base*ones(1,size(TFR_avg_all_in_one,2)));

TFR_left_Z=TFR_avg_all_in_one_Z(1:length(id_lh_dip),index_save);
TFR_right_Z=TFR_avg_all_in_one_Z((length(id_lh_dip)+1):n_grid,index_save);

TFR_left_F=TFR_avg_all_in_one_F(1:length(id_lh_dip),index_save);
TFR_right_F=TFR_avg_all_in_one_F((length(id_lh_dip)+1):n_grid,index_save);
time_labels_save=time_pre_post(index_save);

TFR_sensor_base=mean(TFR_sensor(:,index_base)')';
TFR_sensor_base_std=std(TFR_sensor(:,index_base)')';
TFR_sensor_Z=TFR_sensor-TFR_sensor_base*ones(1,size(TFR_sensor,2));
TFR_sensor_Z=TFR_sensor_Z./(TFR_sensor_base_std*ones(1,size(TFR_sensor,2)));
TFR_sensor_F=TFR_sensor./(TFR_sensor_base*ones(1,size(TFR_sensor,2)));

%********* STEP ************************
% plot the waveforms
fprintf('Plot the result\n')
clf
subplot(2,1,1)
hold on
plot(time_labels_save,TFR_left_F)
plot(time_labels_save,TFR_right_F)
xlabel('Time (ms)')
ylabel('F values')
title('(a) Averaged Time-Freq F waveform')

subplot(2,1,2)
hold on
plot(time_labels_save,TFR_left_Z)
plot(time_labels_save,TFR_right_Z)
xlabel('Time (ms)')
ylabel('Z values')
title('(b) Averaged Time-Freq Z waveform')

drawnow

eval(['print ',result_file,'_F_Z_waveform.jpg -djpeg100'])


eval(['save ',result_file,'.mat']);

fprintf('Write the results to F and Z .stc files...');
write_stc(stc_lh_file_Z,id_lh_dip-1,sf,time_labels_save(1),TFR_left_Z);
write_stc(stc_rh_file_Z,id_rh_dip-1,sf,time_labels_save(1),TFR_right_Z);

write_stc(stc_lh_file_F,id_lh_dip-1,sf,time_labels_save(1),TFR_left_F);
write_stc(stc_rh_file_F,id_rh_dip-1,sf,time_labels_save(1),TFR_right_F);
fprintf('DONE\n')



