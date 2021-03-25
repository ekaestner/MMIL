% cov_from_fif
% Calculate covariance matrix from a fif file
%
% input: fif file
% output: Cvar (covariance matrix)
%
% (c) Mingxiong Huang, Ph.D. 04/10/2007

wd = pwd;			% get current working directory
hfig = figure; 			% establish a new figure window

INPUT = 'cov_from_fif_default_input';	% file to get default inputs
TEMP = 'cov_from_fif_temp_output';	% file to write this run's inputs
set(0,'DefaulttextInterpreter','none');

if(exist([wd '/' INPUT]) == 2), % a defaults file exists from a prior run
  f_default = fopen(INPUT,'r');
else
  f_default = [];		% none found in current directory
end

if(exist([wd '/' TEMP]) == 2),
  delete temp_output		% we'll write new inputs to here
end
f_temp = fopen(TEMP,'w');

datafile=default_s('Name of evoked data set to load (.fif)',f_default,f_temp);
datafile = deblank(datafile);
length_datafile = length(datafile);
stim_num = default_v('Stim Number ',f_default,f_temp);
if strcmp(datafile(length_datafile-3:length_datafile),'.fif') == 1
  fif_data=fiff_read_evoked(datafile,stim_num); 
  sfreq=double(fif_data.info.sfreq); 
  F_orig=double(fif_data.evoked.epochs); 
  t0=double(fif_data.evoked.first)/sfreq; 
  [Coil_306,intpnt_loc_306,intpnt_ori_306]=sensor_306_normal_mne(fif_data);
  F_orig=F_orig(1:length(Coil_306),:);
else
  error('Wrong file name!')
end


flag_project=default_v('Applying SSP? 0 for without SSP, 1 for using SSP: ',f_default,f_temp);

T_project=[];
if flag_project==1,
    temp=fif_data.info.projs;
    for i=1:length(temp) 
        T_project=[T_project temp(i).data.data']; 
    end 
    if size(T_project,1)==102, 
        T_project_tmp=zeros(306,size(T_project,2));
        T_project_tmp((3:3:306),:)=T_project;
        T_project=T_project_tmp;
        clear T_project_tmp
    end    
    T_project=eye(size(T_project,1))-T_project*T_project'; % construct projection <NEW>
else 
    T_project=eye(size(F_orig,1),size(F_orig,1));
end    


id_grad=[];
id_mag=[];
for i=1:length(Coil_306)
    if Coil_306(i).n==2
        id_grad=[id_grad;i];
    elseif Coil_306(i).n==4
        id_mag=[id_mag;i];
    end
end    

time_labels=t0*1000 + [0:(size(F_orig,2)-1)]*1000/sfreq; % in ms

% Now see if user wants to trim this data set
if(size(time_labels) == [1 1]),	% if NOT single time slice then plot
 disp(sprintf('\n You have only one latency: %g ',time_labels));
 begin_lat = time_labels;
 end_lat = time_labels;
else
 F = T_project*F_orig; % plot the sensor waveforms after projection 
 figure(hfig);
 subplot(1,2,1)
 plot(time_labels,F(id_grad,:)'*1e13); 
 title(['Gradiometers in Data set ' datafile]);
 xlabel('Time (ms)');
 ylabel('Response (fT/cm)');
 drawnow
 subplot(1,2,2)
 plot(time_labels,F(id_mag,:)'*1e15); % <NEW>
 title(['magnetometers in Data set ' datafile]);
 xlabel('Time (ms)');
 ylabel('Response (fT)');
 drawnow
 base_lat= default_v('Prestim interval for covariance matrix calculation (ms) [t1 t2] ',f_default,f_temp);
end

% translate to indices
[junk,begin_pre] = min(abs(time_labels - base_lat(1)));
[junk,end_pre] = min(abs(time_labels - base_lat(2)));
base_lat = [begin_pre end_pre];	% adjust user's input for baseline


% now calculate the noise covariance
F_baseline=F(:,base_lat(1):base_lat(2));
mean_base=mean(F_baseline')';
F=F-mean_base*ones(1,size(F,2));   % correct baseline
F_baseline=F(:,base_lat(1):base_lat(2));
Cvar=diag(std(F_baseline').^2);  % diag covariance matrix
navg=1;

cov_filename = default_s('Name of the output file with covariance matrix (.mat) ',f_default,f_temp);
cov_filename = deblank(cov_filename);
length_cov_filename = length(cov_filename);
if strcmp(cov_filename(length_cov_filename-3:length_cov_filename),'.mat') == 1
    % now save the result
    eval(['save ',cov_filename,' Cvar navg']);
else
  error('The output file has to end with .mat!')
end

    
    
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

display(['You result has been saved. Please use xplotter to identify noisy and flat channels!!!'])

