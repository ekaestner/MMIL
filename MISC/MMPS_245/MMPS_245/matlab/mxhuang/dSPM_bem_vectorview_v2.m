% dSPM_bem_vectorview_v2
% dSPM (noise normalized minimum L2-norm) inverse modeling
% for processing VectorView MEG data using BEM forward model
%
% (c) M. Huang 11/08/06
% 
% This program call Dr. Matti Hamalainen's MNE_matlab toolbox
%
% latest change 4/2/2007:
%   modify the lines for loading the covariance matrix from an external
%   file

clear;                                  	% clear working space
INPUT = 'dSPM_vv_bem_project_default_input';	% file to get default inputs
TEMP = 'dSPM_vv_bem_project_temp_output';	% file to write this run's inputs
set(0,'DefaulttextInterpreter','none');
r_cutoff=0.0; % in meter, cutoff the center region 

Figure_window = 'dSPM'; % title of window

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

close all;                      % close all figure windows
hfig = figure; 			% establish a new figure window
set(hfig,'Name',Figure_window);

%**************** STEP ***************************************************
% Step:  Load data into work area.

datafile=default_s('Name of evoked data set to load (.fif)',f_default,f_temp);
datafile = deblank(datafile);
length_datafile = length(datafile);
stim_num = default_v('Stim Number ',f_default,f_temp);
if strcmp(datafile(length_datafile-3:length_datafile),'.fif') == 1
  fif_data=fiff_read_evoked(datafile,stim_num); % load the data structure <NEW>
  sfreq=double(fif_data.info.sfreq); % sampling frequency <NEW>
  F_orig=double(fif_data.evoked.epochs); % the evoked data <NEW>
  t0=double(fif_data.evoked.first)/sfreq; % first time point in sec <NEW>
  [Coil_306,intpnt_loc_306,intpnt_ori_306]=sensor_306_normal_mne(fif_data);
  F_orig=F_orig(1:length(Coil_306),:); % only the MEG data
else
  error('Wrong file name!')
end

% bad channels
bad_chann_names=default_v('List bad MEG channel number, e.g. [2042 1342] ',f_default,f_temp);
bad_chan_id=[];
channames_all=fif_data.info.ch_names; % channel names <NEW>
for i=1:length(channames_all)
    tmp1=channames_all{i};
    for j=1:length(bad_chann_names)
        tmp2=num2str(bad_chann_names(j));
        if length(tmp2)==3
            str_bad1=['MEG 0',tmp2];
            str_bad2=['MEG0',tmp2];
        else
            str_bad1=['MEG ',tmp2];
            str_bad2=['MEG',tmp2];
        end    
        if strcmp(tmp1,str_bad1) | strcmp(tmp1,str_bad2)
            bad_chan_id=[bad_chan_id;i];
        end
    end
end    

for j=1:length(bad_chann_names)
    tmp2=num2str(bad_chann_names(j));
    if length(tmp2)==3
        channames_bad{j}=['MEG0',tmp2];
    else
        channames_bad{j}=['MEG',tmp2];
    end
end    

flag_project=default_v('Applying SSP? 0 for without SSP, 1 for using SSP: ',f_default,f_temp);

T_project=[];
if flag_project==1,
    temp=fif_data.info.projs; % <NEW>
    for i=1:length(temp) % <NEW>
        T_project=[T_project temp(i).data.data']; % <NEW>
    end % <NEW>
    if size(T_project,1)==102, % only magnetometers
        T_project_tmp=zeros(306,size(T_project,2));
        T_project_tmp((3:3:306),:)=T_project;
        T_project=T_project_tmp;
        clear T_project_tmp
    end    
    T_project(bad_chan_id,:)=0.0*T_project(bad_chan_id,:); % ignore the bad channels, <NEW>
    T_project=eye(size(T_project,1))-T_project*T_project'; % construct projection <NEW>
else 
    T_project=eye(size(F_orig,1),size(F_orig,1));
end    


id_grad=[];
id_mag=[];
for i=1:length(Coil_306)
    if isempty(find(bad_chan_id==i))
        if Coil_306(i).n==2
            id_grad=[id_grad;i];
        elseif Coil_306(i).n==4
            id_mag=[id_mag;i];
        end
    end
end    

time_labels=t0*1000 + [0:(size(F_orig,2)-1)]*1000/sfreq; % in ms

% Now see if user wants to trim this data set
if(size(time_labels) == [1 1]),	% if NOT single time slice then plot
 disp(sprintf('\n You have only one latency: %g ',time_labels));
 begin_lat = time_labels;
 end_lat = time_labels;
else
 F = T_project*F_orig; % plot the sensor waveforms after projection <NEW> 
 figure(hfig);
 subplot(1,2,1)
 plot(time_labels,F(id_grad,:)'*1e13); % <NEW>
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
 base_lat= default_v('Prestim interval for baseline (ms) [t1 t2] ',f_default,f_temp);
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
while 1
    cov_method= default_v('Select 1 if Covariance matrix from baseline of this average data set; 2 if from external file: ',f_default,f_temp);
    if cov_method==1,
        Cvar=diag(std(F_baseline').^2);  % diag covariance matrix
        navg=1;
        break;
    elseif cov_method==2,
        cov_mat_name=input('   Input the .mat file name. This file should contains variables Cvar and navg: ','s');
        eval(['load ',cov_mat_name,' Cvar navg']);
        %Cvar=double(noise_matrix.cov_matrix);
        %Cvar=double(noise_cov);
        %navg=input('   The number of trials in the average: ');
        break;
    end
end    


[nlatencies,tmp] = size(time_labels);

subplot(1,2,1)
plot(time_labels,F(id_grad,:)');
title(['Gradiometers in Data set ' datafile]);
xlabel('Time (ms)');
ylabel('Response (fT/cm)');
drawnow
xlim_raw_1=get(gca,'Xlim');
ylim_raw_1=get(gca,'Ylim');

subplot(1,2,2)
plot(time_labels,F(id_mag,:)');
title(['magnetometers in Data set ' datafile]);
xlabel('Time (ms)');
ylabel('Response (fT)');
drawnow
xlim_raw_2=get(gca,'Xlim');
ylim_raw_2=get(gca,'Ylim');


SNR = default_v('SNR of amplitude signal: ',f_default,f_temp);

% now load the dip and dec files
lh_dip_name = default_s('Name of LH Freesurfer .dip file ',f_default,f_temp);
lh_dec_name = default_s('Name of LH Freesurfer .dec file ',f_default,f_temp);
[dip_info_lh,dec_dip_lh]=read_dipdec(lh_dip_name,lh_dec_name);

rh_dip_name = default_s('Name of RH Freesurfer .dip file ',f_default,f_temp);
rh_dec_name = default_s('Name of RH Freesurfer .dec file ',f_default,f_temp);
[dip_info_rh,dec_dip_rh]=read_dipdec(rh_dip_name,rh_dec_name);

% now create the mesh grid
id_lh_dip=find(dec_dip_lh==1);
id_rh_dip=find(dec_dip_rh==1);
grid_mri=[dip_info_lh(1:3,id_lh_dip)';dip_info_rh(1:3,id_rh_dip)']; % in mm
n_grid=size(grid_mri,1);

% now ask the inner skull tri file
tri_file_name = default_s('Name of .tri file for the inner skull surface in MRI coordinate',f_default,f_temp);
[mesh_skull,geo_skull]=load_tri_file(tri_file_name);
temp=max(mesh_skull(:,1))-min(mesh_skull(:,1));
if temp > 100 % original unit in mm
    vertices=mesh_skull(:,1:3)/1000; % mm to m
elseif temp > 10 & temp < 30 % original unit in cm
    vertices=mesh_skull(:,1:3)/100; % cm to m
else
    vertices=mesh_skull(:,1:3); % in m
end
n_vert=size(vertices,1);
faces=geo_skull(:,1:3);


% ask for the .fif file for MRI-to-MEG (HEAD) transformation matrix
mri_fif_name = default_s('Name of .fif file with Tran from MRI to MEG (Head)',f_default,f_temp);
mri_stack=fiff_read_mri(mri_fif_name); % using new program <NEW>
T_head_mri=mri_stack.trans.trans; % the head coor to mri coor <NEW>
T_mri_head=inv(T_head_mri); % the mri coor to head coor <NEW>

% mri to head
grid_head=T_mri_head*[grid_mri'/1000;ones(1,n_grid)]; % in meters
grid_head=grid_head(1:3,:)'; %  in meter

vertices_head=T_mri_head*[vertices';ones(1,n_vert)];
vertices_head=vertices_head(1:3,:)'; 

% now fit the center of sphere
cen_sph = fminsearch(@(cen_sph) sphererr(cen_sph,vertices_head),mean(vertices_head)')';

% shift the grid
grid_head_cen=grid_head-ones(n_grid,1)*cen_sph;

% remove the ones that are too close to the center
grid_head_cen_rtp=cart2rtp(grid_head_cen);
id_keep=find(grid_head_cen_rtp(:,1)>r_cutoff); % find the points that are NOT too close to the center (20 mm)
grid_head_cen_keep=grid_head_cen(id_keep,:);
grid_head_keep=grid_head(id_keep,:);

% ask the LU file name
H_inv_filename = default_s('File name (.mat) of pre-calculated BEM  (will create one if not exist)',f_default,f_temp);

% now draw the BEM mesh, source nodes, and sensor locations
clf
bem_mesh_plot({vertices_head*1000},{faces})
plot3d(grid_head*1000,'b*');
plot3d(intpnt_loc_306*1000,'go')
title(['Sensor and source grid ' datafile]);
xlabel('x (mm)'); 
ylabel('y (mm)'); 
zlabel('z (mm)'); 
hold off
drawnow

fprintf('Now, check your dipole grid, bem mesh, and sensor location plot\n');

% pick up the gradiometer, magnetometer, or both
analyze_mode = default_v('Select 1 for gradiometers, 2 for magnetometers, 3 for both ',f_default,f_temp);


% ask for summary file name
summary_file = default_s('Name of output summary file ',f_default,f_temp);


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


% put together for BEM
bem_input.R_meg=intpnt_loc_306; 
bem_input.O_meg=intpnt_ori_306;
bem_input.R_eeg=[];
bem_input.vertices={vertices_head};
bem_input.faces={faces};
bem_input.sigma=0.3; % conductivity
bem_input.mode=2; % MEG only
bem_input.basis_opt=1; % linear potential function
bem_input.test_opt=0; % collocation
bem_input.ISA=0; % Inhibit Isolated Skull Approach
bem_input.fn_eeg='bem_eeg_temp';
bem_input.fn_meg='bem_meg_temp';

% calculate transfer matrix using BEM
%return_var = bem_transfer_1shell(bem_input,[],1);
return_var = bem_transfer_3shell(bem_input,[],H_inv_filename,1);
G_bem_xyz_inp=bem_gainmat(grid_head_keep,return_var,0);
G_bem_xyz=T_project*gain_intpnt2chan(Coil_306,G_bem_xyz_inp.meg);

for i=1:length(Coil_306)
    if Coil_306(i).n==2
        G_bem_xyz(i,:)=G_bem_xyz(i,:)*1e-13; % into T/m
    elseif Coil_306(i).n==4
        G_bem_xyz(i,:)=G_bem_xyz(i,:)*1e-15; % into T
    end
end

clear G_bem_xyz_inp
%G_bem_tp=gain_xyz_tp(G_bem_xyz,grid_head_cen_keep);

% now obtain the two strong orientations
%G_bem_b2=zeros(size(G_bem_xyz,1),2*length(id_keep));
%ori_b2=zeros(3,2*length(id_keep));
%for i=1:length(id_keep)
%    id_temp1=[-1:0]+2*i;
%    id_temp2=[-2:0]+3*i;
%    [U_temp,S_temp,V_temp]=svd(G_bem_xyz(:,id_temp2));
%    G_bem_b2(:,id_temp1)=G_bem_xyz(:,id_temp2)*V_temp(:,1:2);
%    ori_b2(:,id_temp1)=V_temp(:,1:2);
%end    
%clear G_bem_xyz

%*************** STEP ****************************************************
% Step: now start the dSPM norm calculation

fprintf('Start the dSPM calculation...');

if analyze_mode==1 % only fit the gradiometer
    id_select=id_grad;
elseif analyze_mode==2 % only fir the magnetometer
    id_select=id_mag;
else
    id_select=[id_grad;id_mag];
end  

data_input.meg={F(id_select,:)}; % cell array
data_input.t=time_labels;
data_input.ncov=Cvar(id_select,id_select);
data_input.navg=navg;
data_input.SNR=SNR;
R=sparse(3*length(id_keep),3*length(id_keep));
for i=1:3*length(id_keep),
    R(i,i)=1;
end
data_input.R=R;

[dSPM_F_keep,J_dip_moment_keep,junk,mn_inverse]=mn_inverse_operator(data_input,G_bem_xyz(id_select,:));

fprintf('DONE\n');

% the predicted fields
F_predict=G_bem_xyz*J_dip_moment_keep{1};      

sqrt_dSPM_F=zeros(n_grid,length(time_labels));
sqrt_dSPM_F(id_keep,:)=sqrt(dSPM_F_keep{1});

% now write the .stc file
if exist([wd,'/stc'])~=7
    !mkdir stc
end    
fprintf('Write the dSPM results to .stc files...');
stc_lh_file=['stc/',summary_file,'-lh.stc'];
stc_rh_file=['stc/',summary_file,'-rh.stc'];

sqrt_dSPM_F_left=sqrt_dSPM_F(1:length(id_lh_dip),:);
sqrt_dSPM_F_right=sqrt_dSPM_F((length(id_lh_dip)+1):n_grid,:);

write_stc(stc_lh_file,id_lh_dip-1,sfreq,time_labels(1),sqrt_dSPM_F_left);
write_stc(stc_rh_file,id_rh_dip-1,sfreq,time_labels(1),sqrt_dSPM_F_right);
fprintf('DONE\n')

%*************** STEP **********************************************
% now sort the index sources according to the F values
norm_sqrt_dSPM_F=rownorm(sqrt_dSPM_F);
[dum,id_F_sort]=sort(norm_sqrt_dSPM_F);
id_F_sort=flipud(id_F_sort);

np=size(sqrt_dSPM_F,1);
n_left=length(id_lh_dip);
n_right=length(id_rh_dip);

id_fsurfer_vertex=ones(np,2);

for i=1:np
   if id_F_sort(i)<=n_left
      id_fsurfer_vertex(i,:)=[id_lh_dip(id_F_sort(i))-1 1];
   else
      id_fsurfer_vertex(i,:)=[id_rh_dip(id_F_sort(i)-n_left)-1 2];
   end
end  

%********* STEP ************************
% plot the waveforms

fprintf('Plot the result\n')
clf
subplot(2,2,1)
plot(time_labels,F(id_grad,:))
xlabel('Time (ms)')
ylabel('B (fT/cm)')
title('(a) Measured gradiometer waveform')
Xlim=get(gca,'Xlim');
Ylim=get(gca,'Ylim');

subplot(2,2,2)
plot(time_labels,F_predict(id_grad,:))
xlabel('Time (ms)')
ylabel('B (fT/cm)')
title('(b) Predicted gradiometer waveform')
set(gca,'Xlim',Xlim,'Ylim',Ylim)

subplot(2,2,3)
plot(time_labels,F(id_grad,:)-F_predict(id_grad,:))
xlabel('Time (ms)')
ylabel('B (fT/cm)')
title('(c) Residual gradiometer waveform')
set(gca,'Xlim',Xlim,'Ylim',Ylim)

subplot(2,2,4)
plot(time_labels,sqrt_dSPM_F)
xlabel('Time (ms)')
ylabel('sqrt of F values')
title('(d) sqrt of F waveform')

drawnow

figure(2)
clf
subplot(2,2,1)
plot(time_labels,F(id_mag,:))
xlabel('Time (ms)')
ylabel('B (fT)')
title('(a) Measured magnetometer waveform')
Xlim=get(gca,'Xlim');
Ylim=get(gca,'Ylim');

subplot(2,2,2)
plot(time_labels,F_predict(id_mag,:))
xlabel('Time (ms)')
ylabel('B (fT)')
title('(b) Predicted magnetometer waveform')
set(gca,'Xlim',Xlim,'Ylim',Ylim)

subplot(2,2,3)
plot(time_labels,F(id_mag,:)-F_predict(id_mag,:))
xlabel('Time (ms)')
ylabel('B (fT)')
title('(c) Residual magnetometer waveform')
set(gca,'Xlim',Xlim,'Ylim',Ylim)

subplot(2,2,4)
plot(time_labels,sqrt_dSPM_F)
xlabel('Time (ms)')
ylabel('sqrt of F values')
title('(d) sqrt of F waveform')


drawnow


clear return_var G_bem_xyz dSPM_F_keep

fprintf('Save result and figures...')

eval(['save ',summary_file,'.mat'])

figure(1)
eval(['print ',summary_file,'_grad_waveform.jpg -djpeg100'])
figure(2)
eval(['print ',summary_file,'_mag_waveform.jpg -djpeg100'])

fprintf('DONE\n');
fprintf('Thank you!\n') 
