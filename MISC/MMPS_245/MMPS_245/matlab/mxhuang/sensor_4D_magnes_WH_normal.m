function [Coil,intpnt_loc,intpnt_ori]=sensor_4D_magnes_WH_normal(fif_filename,apos_filename);

% function [Coil,intpnt_loc,int_pnt_ori]=sensor_4D_magnes_WH_normal(fif_filename,apos_filename);
% calling without fif_filename will return device coordinates
%
% SENSOR_4D_MAGNES_WH_NORMAL Return the Vector-view 306 
% coil information in NORMAL FORMAT
%
% Coil: 1 by 248 structure. 
%       Coil(i).n the number of integration points in each channel
%       Coil(i).wei are the weight matrix for the all integration points in channel i; 
%       Coil(i).loc location of integration points;
%       Coil(i).ori orientation of the integration
% intpnt_loc: location of all integration points concatenated 
% intpnt_ori: orientation of all integration points concatenated

if(nargin==0),
    load sensor_248_default.mat
    T_device_head=eye(4);
elseif exist(fif_filename)==0
    error('Fif file does not exist')
else
    [sensor.name,sensor.kind,sensor.logname,sensor.type,sensor.loc]=channames(fif_filename);
    T_device_head=load_trans_apos(apos_filename,0);
    %T_device_head=eye(4);
end

Coil=[];
intpnt_loc=[];
intpnt_ori=[];


for i=1:length(sensor.name),
   Trans=sensor.loc{i};
   if strcmp(sensor.name{i}(1:3),'MEG')==1,     % now MEG coils
       Rcenter=Trans(1:3,4)';
       Ux=Trans(1:3,1)';
       Uy=Trans(1:3,2)';
       Uz=Trans(1:3,3)';
       Coil_ori_tmp=[Uz; Uz];	% directions
       if sensor.type(i)==4001, 
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
