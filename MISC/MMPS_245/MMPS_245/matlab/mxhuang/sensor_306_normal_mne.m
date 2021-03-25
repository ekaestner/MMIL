function [Coil,intpnt_loc,intpnt_ori]=sensor_306_normal_mne(fif_data);

% function [Coil,intpnt_loc,int_pnt_ori]=sensor_306_normal_mne(fif_filename);
%
% SENSOR_306_NORMAL Return the Vector-view 306 
% coil information in NORMAL FORMAT
%
% fif_data is the data structure returned from 'fiff_read_evoked.m' from
% Dr. Matti Hamalainen's MNE_matlab toolbox
%
% Coil: 1 by 306 structure. 
%       Coil(i).n the number of integration points in each channel
%       Coil(i).wei are the weight matrix for the all integration points in channel i; 
%       Coil(i).loc location of integration points;
%       Coil(i).ori orientation of the integration
% intpnt_loc: location of all integration points concatenated 
% intpnt_ori: orientation of all integration points concatenated
%
% (c) Mingxiong Huang, Ph.D.
% last update 11/08/2006

Coil=[];
intpnt_loc=[];
intpnt_ori=[];

T_device_head=fif_data.info.dev_head_t.trans;

for i=1:length(fif_data.info.ch_names),
   Trans=reshape(fif_data.info.chs(i).loc,3,4);
   Trans=[Trans(:,2:4) Trans(:,1);0 0 0 1];
   if strcmp(fif_data.info.ch_names{i}(1:3),'MEG')==1,     % now MEG coils
       Rcenter=Trans(1:3,4)';
       Ux=Trans(1:3,1)';
       Uy=Trans(1:3,2)';
       Uz=Trans(1:3,3)';
       if fif_data.info.chs(i).coil_type==3012 | fif_data.info.chs(i).coil_type==3013,
           Coil(i).n=2;
           Coil(i).wei=diag([1 -1]/1.68); % fT/cm
           Coil_loc_tmp=[(Rcenter+[8.4 0 0.3]*[Ux;Uy;Uz]/1000);...
                        (Rcenter+[-8.4 0 0.3]*[Ux;Uy;Uz]/1000)];
           Coil_loc_tmp=T_device_head*[Coil_loc_tmp';[1 1]];
           Coil(i).loc=Coil_loc_tmp(1:3,:)';
           Coil_ori_tmp=[Uz; Uz];	% directions
           Coil_ori_tmp=T_device_head(1:3,1:3)*Coil_ori_tmp';
           Coil(i).ori=Coil_ori_tmp(1:3,:)';
           intpnt_loc=[intpnt_loc;Coil(i).loc];
           intpnt_ori=[intpnt_ori;Coil(i).ori];
       elseif fif_data.info.chs(i).coil_type==3022 | fif_data.info.chs(i).coil_type==3023,
           Coil(i).n=4;
           Coil(i).wei=diag([1 1 1 1]/4); % fT
           Coil_loc_tmp=[(Rcenter+[6.45 6.45 0.3]*[Ux;Uy;Uz]/1000);... 
                        (Rcenter+[6.45 -6.45 0.3]*[Ux;Uy;Uz]/1000);...
                        (Rcenter+[-6.45 6.45 0.3]*[Ux;Uy;Uz]/1000);... 
                        (Rcenter+[-6.45 -6.45 0.3]*[Ux;Uy;Uz]/1000);];
           Coil_loc_tmp=T_device_head*[Coil_loc_tmp';[1 1 1 1]];
           Coil(i).loc=Coil_loc_tmp(1:3,:)';
           Coil_ori_tmp=[Uz;Uz;Uz;Uz];	% directions
           Coil_ori_tmp=T_device_head(1:3,1:3)*Coil_ori_tmp';
           Coil(i).ori=Coil_ori_tmp(1:3,:)';
           intpnt_loc=[intpnt_loc;Coil(i).loc];
           intpnt_ori=[intpnt_ori;Coil(i).ori];      
       elseif fif_data.info.chs(i).coil_type==3024, 
           Coil(i).n=4;
           Coil(i).wei=diag([1 1 1 1]/4); % fT
           Coil_loc_tmp=[(Rcenter+[5.25 5.25 0.3]*[Ux;Uy;Uz]/1000);... 
                        (Rcenter+[5.25 -5.25 0.3]*[Ux;Uy;Uz]/1000);...
                        (Rcenter+[-5.25 5.25 0.3]*[Ux;Uy;Uz]/1000);... 
                        (Rcenter+[-5.25 -5.25 0.3]*[Ux;Uy;Uz]/1000);];
           Coil_loc_tmp=T_device_head*[Coil_loc_tmp';[1 1 1 1]];
           Coil(i).loc=Coil_loc_tmp(1:3,:)';
           Coil_ori_tmp=[Uz;Uz;Uz;Uz];	% directions 
           Coil_ori_tmp=T_device_head(1:3,1:3)*Coil_ori_tmp';
           Coil(i).ori=Coil_ori_tmp(1:3,:)';     
           intpnt_loc=[intpnt_loc;Coil(i).loc];
           intpnt_ori=[intpnt_ori;Coil(i).ori]      
       end
   end
end   
     

return
