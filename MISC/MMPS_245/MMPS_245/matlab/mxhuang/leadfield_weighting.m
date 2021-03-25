% leadfield_weighting
% This program calculate the weighting factor for points at different 
% head tissue boundary, with respect to the VectorVirew MEG sensor
% configuration.
% It requres the MRI coordinates for the following landmarks
% 1) Left PA
% 2) Nasion
% 3) Right PA
% and of course the MRI coordinated of the points of interested
%
% all these coordinates can be obtained from Freesurfer
% (c) Mingxiong Huang, Ph.D.

clear;                                  	% clear working space
INPUT = 'leadfield_weighting_default_input';	% file to get default inputs
TEMP = 'leadfield_weighting_temp_output';	% file to write this run's inputs
set(0,'DefaulttextInterpreter','none');

Figure_window = 'Leadfield Weighting'; % title of window

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


% now getting the MRI coordinates 
r_LPA=default_v('MRI coordinates for Left PA (mm) ',f_default,f_temp);
r_NAS=default_v('MRI coordinates for Nasion (mm) ',f_default,f_temp);
r_RPA=default_v('MRI coordinates for Rigth PA (mm) ',f_default,f_temp);

% now the MRI coordinated for points of interest
fprintf('----------------------------------------------\n'); 
r=default_m('MRI coordinates for points of interest (mm) ',f_default,f_temp);

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


nr=size(r,1);
w=zeros(nr,1);
r_cen=r-ones(nr,1)*[0 0 50]; % shift to the center of helmet

% now load the MEG sensor configuration
[Coil,intpnt_loc,intpnt_ori]=sensor_306_normal_device;
n_coil=size(intpnt_loc,1);
plot3d(intpnt_loc*1000,'ro')
hold on
plot3d(r_cen,'*')
hold off

fprintf('\n\tx(mm)\ty(mm)\tz(mm)\tweight\n')
for i=1:nr
    r_dip=r_cen(i,:)/1000; % mm to meter
    r_minum_rp=intpnt_loc-ones(n_coil,1)*r_dip;
    w(i)=mean(ones(n_coil,1)./(rownorm(r_minum_rp).^2));
    fprintf('\t%.f\t%.f\t%.f\t%.f\n',r(i,1),r(i,2),r(i,3),w(i));
end    