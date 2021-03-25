function [Coil,intpnt_loc,intpnt_ori]=sensor_4D_magnes_WH_normal_mne(fif_data);

% function [Coil,intpnt_loc,int_pnt_ori]=sensor_4D_magnes_WH_normal_mne(fif_data);
% calling without fif_filename will return device coordinates
%
% SENSOR_4D_MAGNES_WH_NORMAL Return the 4D NeuroImaging Magnes 
% coil information in NORMAL FORMAT
%
% fif_data is the data structure returned from 'fiff_setup_read_raw.m' from
% Dr. Matti Hamalainen's MNE_matlab toolbox
%
%
% Coil: 1 by 248 structure. 
%       Coil(i).n the number of integration points in each channel
%       Coil(i).wei are the weight matrix for the all integration points in channel i; 
%       Coil(i).loc location of integration points;
%       Coil(i).ori orientation of the integration
% intpnt_loc: location of all integration points concatenated 
% intpnt_ori: orientation of all integration points concatenated
% 
% (c) Mingxiong Huang, Ph.D.
% Last changes: 11/12/2006


Coil=[];
intpnt_loc=[];
intpnt_ori=[];

T_device_head=fif_data.info.dev_head_t.trans;


for i=1:length(fif_data.info.ch_names),
   Trans=reshape(fif_data.info.chs(i).loc,3,4);
   Trans=[Trans(:,2:4) Trans(:,1);0 0 0 1];
   if strcmp(fif_data.info.ch_names{i}(1:3),'MEG')==1, % now MEG coils
       Rcenter=Trans(1:3,4)';
       Ux=Trans(1:3,1)';
       Uy=Trans(1:3,2)';
       Uz=Trans(1:3,3)';
       Coil_ori_tmp=[Uz; Uz];	% directions
       if fif_data.info.chs(i).coil_type==4001, 
           Coil(i).n=4;
           Coil(i).wei=diag([1 1 1 1]/4); % fT
           Coil_loc_tmp=[(Rcenter+[5.75 5.75 0.0]*[Ux;Uy;Uz]/1000);... 
                        (Rcenter+[5.75 -5.75 0.0]*[Ux;Uy;Uz]/1000);...
                        (Rcenter+[-5.75 5.75 0.0]*[Ux;Uy;Uz]/1000);... 
                        (Rcenter+[-5.75 -5.75 0.0]*[Ux;Uy;Uz]/1000);];
           Coil_loc_tmp=T_device_head*[Coil_loc_tmp';[1 1 1 1]];
           Coil(i).loc=Coil_loc_tmp(1:3,:)';
           Coil_ori_tmp=[Uz;Uz;Uz;Uz];	% directions 
           Coil_ori_tmp=T_device_head(1:3,1:3)*Coil_ori_tmp';
           Coil(i).ori=Coil_ori_tmp(1:3,:)';     
           intpnt_loc=[intpnt_loc;Coil(i).loc];
           intpnt_ori=[intpnt_ori;Coil(i).ori];      
       end
   end
end   
     

return
