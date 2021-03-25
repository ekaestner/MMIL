function MEG_info = ts_prep_MEG_info(avg_data)
%function MEG_info = ts_prep_MEG_info(avg_data)
%
% Required Input:
%   avg_data: structure containing the following fields:
%      num_sensors    (int)
%      sfreq          (double)
%      sensor_info    (struct)
%        type         (string)
%        label        (string)
%        loc          (4x4 matrix of doubles)
%        badchan      (int: 0 or 1)
%      coor_trans     (struct)
%        device2head  (4x4 matrix of doubles) (or empty)
%        mri2head     (4x4 matrix of doubles) (or empty)
%      averages       (struct)
%        event_code   (int)
%        num_trials   (int)
%        num_rejects  (struct)
%          mag        (int)
%          grad       (int)
%          eeg        (int)
%          eog        (int)
%          manual     (int)
%        time         (vector of doubles) (in msec)
%        data         (matrix of doubles) (num_sensors x num samples)
%        stdev        (matrix of doubles) (num_sensors x num samples)
%      noise          (struct)
%        num_trials   (int)
%        scale_fact   (double)
%        covar        (matrix of doubles) (num_sensors x num_sensors)
%
% Note: minimally necessary information is "sensor_info" and "coor_trans"
%     fields of avg_data structure, which is returned by ts_avg_fif_data
%
% Output:
%   MEG 306 coil information
%   Note: coordinates are returned in "head" space
%
%   MEG_info.Coil: 1 by 306 struct array
%       Coil(i).n the number of integration points in each channel
%       Coil(i).wei are the weight matrix for the all integration points in channel i; 
%       Coil(i).loc location of integration points;
%       Coil(i).ori orientation of the integration
%   MEG_info.intpnt_loc: location of all integration points concatenated 
%   MEG_info.intpnt_ori: orientation of all integration points concatenated
%
%
% Based on sensor_306_normal By M.X. Huang, PhD, March 2004
%
% Created:  07/31/06 by Don Hagler
% Last Mod: 10/11/13 by Don Hagler
%

MEG_info=[];
Coil=[];
intpnt_loc=[];
intpnt_ori=[];
labels={};

T_device_head=avg_data.coor_trans.device2head;

for i=1:avg_data.num_sensors;
  loc = avg_data.sensor_info(i).loc;
  if strcmp(upper(avg_data.sensor_info(i).label(1:3)),'MEG')==1,     % now MEG coils
    Rcenter=loc(1:3,4)';
    Ux=loc(1:3,1)';
    Uy=loc(1:3,2)';
    Uz=loc(1:3,3)';
    if avg_data.sensor_info(i).type==3012 | avg_data.sensor_info(i).type==3013,
        Coil(i).n=2;
        Coil(i).wei=diag([1 -1]/1.68); % fT/cm
        Coil_loc_tmp=[(Rcenter+[8.4 0 0.3]*[Ux;Uy;Uz]/1000);...
                     (Rcenter+[-8.4 0 0.3]*[Ux;Uy;Uz]/1000)];
        Coil_loc_tmp=T_device_head*[Coil_loc_tmp';[1 1]];
        Coil_ori_tmp=[Uz; Uz];	% directions
    elseif avg_data.sensor_info(i).type==6001, % axial grad on KIT; added by Andrei Irimia, 06/25/09
        Coil(i).n=2;
        Coil(i).wei=diag([1 1]/5); % fT/cm
        Coil_loc_tmp=[(Rcenter+[7.75 0 0]*[Ux;Uy;Uz]/1000);...
                      (Rcenter+[7.75 0 0]*[Ux;Uy;Uz]/1000)];                
        Coil_loc_tmp=T_device_head*[Coil_loc_tmp';[1 1]];
        Coil_ori_tmp=[Uz; Uz];	% directions
    elseif avg_data.sensor_info(i).type==3022 | avg_data.sensor_info(i).type==3023,
        Coil(i).n=4;
        Coil(i).wei=diag([1 1 1 1]/4); % fT
        Coil_loc_tmp=[(Rcenter+[6.45 6.45 0.3]*[Ux;Uy;Uz]/1000);... 
                     (Rcenter+[6.45 -6.45 0.3]*[Ux;Uy;Uz]/1000);...
                     (Rcenter+[-6.45 6.45 0.3]*[Ux;Uy;Uz]/1000);... 
                     (Rcenter+[-6.45 -6.45 0.3]*[Ux;Uy;Uz]/1000);];
        Coil_loc_tmp=T_device_head*[Coil_loc_tmp';[1 1 1 1]];
        Coil_ori_tmp=[Uz;Uz;Uz;Uz];	% directions
    elseif avg_data.sensor_info(i).type==3024, 
        Coil(i).n=4;
        Coil(i).wei=diag([1 1 1 1]/4); % fT
        Coil_loc_tmp=[(Rcenter+[5.25 5.25 0.3]*[Ux;Uy;Uz]/1000);... 
                     (Rcenter+[5.25 -5.25 0.3]*[Ux;Uy;Uz]/1000);...
                     (Rcenter+[-5.25 5.25 0.3]*[Ux;Uy;Uz]/1000);... 
                     (Rcenter+[-5.25 -5.25 0.3]*[Ux;Uy;Uz]/1000);];
        Coil_loc_tmp=T_device_head*[Coil_loc_tmp';[1 1 1 1]];
        Coil_ori_tmp=[Uz;Uz;Uz;Uz];	% directions 
    elseif avg_data.sensor_info(i).type==4001, 
        Coil(i).n=4;
        Coil(i).wei=diag([1 1 1 1]/4); % fT
        Coil_loc_tmp=[(Rcenter+[5.75 5.75 0.0]*[Ux;Uy;Uz]/1000);... 
                     (Rcenter+[5.75 -5.75 0.0]*[Ux;Uy;Uz]/1000);...
                     (Rcenter+[-5.75 5.75 0.0]*[Ux;Uy;Uz]/1000);... 
                     (Rcenter+[-5.75 -5.75 0.0]*[Ux;Uy;Uz]/1000);];
        Coil_loc_tmp=T_device_head*[Coil_loc_tmp';[1 1 1 1]];
        Coil_ori_tmp=[Uz;Uz;Uz;Uz];	% directions 
    end
    Coil(i).loc=Coil_loc_tmp(1:3,:)';
    Coil_ori_tmp=T_device_head(1:3,1:3)*Coil_ori_tmp';
    Coil(i).ori=Coil_ori_tmp(1:3,:)';     
    intpnt_loc=[intpnt_loc;Coil(i).loc];
    intpnt_ori=[intpnt_ori;Coil(i).ori];      
    labels{i} = avg_data.sensor_info(i).label;
  end
end

MEG_info.Coil = Coil;
MEG_info.intpnt_loc = intpnt_loc;
MEG_info.intpnt_ori = intpnt_ori;
MEG_info.labels = labels;

return

