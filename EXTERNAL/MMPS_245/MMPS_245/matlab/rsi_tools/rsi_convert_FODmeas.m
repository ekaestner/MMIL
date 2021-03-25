function rsi_convert_FODmeas(FODmeas,fstem,varargin)
%function rsi_convert_FODmeas(FODmeas,fstem,[options])
%
% Usage:
%  rsi_convert_FODmeas(MFmeas,fstem,'key1', value1,...);
%
% Required Parameters:
%   MFmeas: structure containing multi-FOD measures
%     volmask : dilated brain mask         size = [nx,ny,nz]  
%     volsz   : size of original data           = [nx,ny,nz,nf]
%     volT    :  norm of all parameters    size = [nx,ny,nz]
%     volF0   : 0th order FOD              size = [nx,ny,nz]
%     volN0   : volF0 normalized by volT   size = [nx,ny,nz]
%     volF2   : norm of 2nd order FOD      size = [nx,ny,nz]
%     volN2   : volF2 normalized by volT   size = [nx,ny,nz]
%     volF4   : norm of 4th order FOD      size = [nx,ny,nz]
%     volN4   : volF4 normalized by volT   size = [nx,ny,nz]
%     volFD   : norm of directional FOD    size = [nx,ny,nz,ns]
%               (F2 and F4 combined)
%     volND   : volFD normalized by volT   size = [nx,ny,nz,ns]
%     volFT   : norm of total FOD          size = [nx,ny,nz,ns]
%               (F0, F2, and F4 combined)
%     volNT   : volFT normalized by volT   size = [nx,ny,nz,ns]
%     volV0   : orientation of maximum FOD size = [nx,ny,nz,3]
%     volAU   : angular uncertainty        size = [nx,ny,nz]
%   fstem: output file stem (include full path)
%
% Optional Input Parameters:
%   'regT1flag': [0|1|2] whether to resample DT calculations into T1-alignment
%     This is done in addition to converting in DTI resolution
%     0: do not resample DT results to T1 resolution
%     1: apply DTI to T1 registration but keep DTI resolution,
%        (padded with extra slices to match number of in-plane voxels)
%     2: apply DTI to T1 registration and resample to T1 resolution
%     {default = 0}
%   'zpadflag': [0|1] pad slices to match number of in-plane voxels
%     avoids clipping of image after registration to T1
%     (only applicable if regT1flag = 1)
%     {default = 1}
%   'RegInfo': structure with fields (Required if regT1flag > 0):
%     M_T1_to_T2 : transformation matrix
%     M_T2 : DTI vox2ras matrix
%     M_T1 : T1 vox2ras matrix
%     volsz_T2 : vector of DTI voxel sizes [x,y,z]
%     volsz_T1 : vector of T1 voxel sizes [x,y,z]
%     {default = []}
%   'M': vox2ras matrix of input diffusion data
%     ignored if RegInfo is supplied
%     {default = eye(4)}
%   'measlist': cell array of DT measures to convert to mgh/mgz
%     {default = {'T','F0','N0','F2','N2','F4','N4','FD','ND',...
%                 'FT','NT','V0','AU'}}
%   'mgz_flag': [0|1] use compressed mgz format (otherwise mgh)
%     {default = 1}
%   'dat_flag': convert DT measures to dat format (for DTI Studio)
%     {default = 0}
%   'dat_measlist': cell array of DT measures to convert to dat
%     Must be included in 'measlist'
%     {default = {'F2','V0'}
%   'dat_frame': frame number of multi-frame volumes to convert to dat
%     Applies to volumes like F0, F2, F4, not V0
%     {default = 1}
%   'nii_flag': convert DT measures to nii format (for FSL)
%     {default = 0}
%   'nii_measlist': cell array of DT measures to convert to nii
%     must be included in 'measlist'
%     {default = {'T','F0','N0','F2','N2','F4','N4','FD','ND',...
%                 'FT','NT','V0','AU'}}
%   'nii_out_orient': output orientation for nii files
%     e.g. 'RAS', 'LPI', etc.
%     If empty, keep original orientation
%     {default = []}
%   'verbose_flag': [0|1] display status messages
%     {default = 0}
%   'forceflag': run calculations even if output files exist
%     {default = 0}
%
% Created:  05/18/12 by Don Hagler
% Last Mod: 07/31/15 by Don Hagler
%

%% todo: (F2-F4)/(F2+F4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'regT1flag',0,[0:2],...
  'zpadflag',true,[false true],...
  'M',eye(4),[],...
  'RegInfo',[],[],...
  'measlist',{'T','F0','N0','F2','N2','F4','N4','FD','ND',...
              'FT','NT','V0','AU'},[],...
  'mgz_flag',true,[false true],...
  'dat_flag',false,[false true],...
  'dat_measlist',{'F2','V0'},[],...
  'dat_frame',1,[],...
  'nii_flag',false,[false true],...
  'nii_measlist',{'T','F0','N0','F2','N2','F4','N4','FD','ND',...
                  'FT','NT','V0','AU'},[],...
  'nii_out_orient',[],[],...
  'verbose_flag',false,[false true],...
  'forceflag',false,[false true],...
...
  'dat_revsliceflag',true,[false true],...
  'dat_permvec',[1,2,3],[],...  
});

args = mmil_parms2args(parms);
rsi_convert_MFmeas(FODmeas,fstem,args{:})


