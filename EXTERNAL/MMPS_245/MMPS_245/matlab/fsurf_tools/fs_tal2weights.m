function [fnames_weights,fnames_roinames] = fs_tal2weights(fname_points,varargin)
%function [fnames_weights,fnames_roinames] = fs_tal2weights(fname_points,[options])
%
% Purpose: sample Talairach coordinates in a csv file
%   to a subject's cortical surfaces and create a weighted ROI file
%
% Required Input:
%   fname_points: csv file containing Talairach coordinates
%     required column headers: 'x', 'y', 'z'
%     optional column headers:
%       'roiname': may choose different column with 'namelabel' option
%           if not supplied, ROIs will be given numerically sequential names
%       'hemi': must be 'lh' or 'rh'
%           if not supplied, will infer from sign of x
%             if x >= -hemi_thresh, then right
%             if x <= +hemi_thresh, then left
%       'MNIflag': may choose different column name with 'MNIlabel' option
%           overrides 'MNIflag' option for individual points
%     additional columns will be ignored
%     
% Optional Input:
%   'outdir': output directory
%     {default = pwd}
%   'outstem': output file stem
%     {default = 'weights'}
%   'subj': subject name
%     {default = 'fsaverage'}
%   'subjdir': FreeSurfer subject root directory
%     {default = $FREESURFER_HOME/subjects}
%   'namelabel': header of column containing ROI names
%     {default = 'roiname'}
%   'MNIflag': [0|1] whether coordinates are MNI or "true" Talairach
%     {default = 1}
%   'MNIlabel': header of column containing MNI flags for each point
%     overrides 'MNIflag' option
%     {default = 'MNIflag'}
%   'hemi_thresh': threshold for assigning hemisphere
%     values greater than 0 allow mid-line points to be assigned to both
%     {default = 0}
%   'vol_sm': smoothing sigma in volume
%     {default = 0}
%   'surf_sm': number of smoothing iterations on surface
%     {default = 50}
%   'thresh': threshold applied after smoothing
%     {default = 0.1}
%   'binarize_flag': create binarized mask instead of continuous values
%     {default = 0}
%   'hemilist': cell array of cortical hemispheres to search
%     {default = {'lh','rh'}}
%   'verbose': [0|1] display status messages
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%     
% Output:
%   fnames_weights: cell array of multi-frame weights files created
%   fnames_roinames: cell array of text files containing ROI names
%
% Created:  03/03/13 by Don Hagler
% Last Mod: 11/13/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'outdir',pwd,[],...
  'outstem','weights',[],...
  'subj','fsaverage',[],...
  'subjdir',[],[],...
  'namelabel','roiname',[],...
  'MNIflag',true,[false true],...
  'MNIlabel','MNIflag',[],...
  'hemi_thresh',0,[0,10],...
  'vol_sm',0,[0,1000],...
  'surf_sm',50,[0,1000],...
  'thresh',0.1,[0,1],...
  'binarize_flag',false,[false true],...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'verbose',true,[false true],...
  'forceflag',false,[false true],...
...
  'required_fields',{'x','y','z'},[],...
  'hemilabel','hemi',[],...
  'tmpdir',[],[],...
  'outext','.mgz',{'.mgh','.mgz'},...
...
  'surf_tags',{'subj','subjdir','MNIflag','surfname','vol_sm','surf_sm',...
               'thresh','binarize_flag','verbose','tmpdir'},[],...
});
fnames = {};

if isempty(parms.tmpdir)
  parms.tmpdir = [parms.outdir '/tmp'];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sample points to surface for each hemisphere
for h=1:length(parms.hemilist)
  hemi = parms.hemilist{h};
  fname_out = sprintf('%s/%s-%s%s',...
    parms.outdir,parms.outstem,hemi,parms.outext);
  fname_csv = sprintf('%s/%s-%s-info.csv',...
    parms.outdir,parms.outstem,hemi);
  fname_txt = sprintf('%s/%s-%s-roinames.txt',...
    parms.outdir,parms.outstem,hemi);
  if ~exist(fname_out,'file') ||...
     ~exist(fname_csv,'file') ||...
     ~exist(fname_txt,'file') || parms.forceflag

    if parms.verbose
      fprintf('%s: loading points file %s...\n',...
        mfilename,fname_points);
    end;
    roi_info = mmil_csv2struct(fname_points);
    fstrs = fieldnames(roi_info);

    % check for required fields
    if length(intersect(parms.required_fields,fstrs))<length(parms.required_fields)
      error('missing one or more required fields ( %s) in %s',...
        sprintf('%s ',parms.required_fields),fname_points);
    end;

    % exclude any points with missing coordinates
    xtmp = {roi_info.x};
    ytmp = {roi_info.y};
    ztmp = {roi_info.z};
    ind_excl = find(cellfun(@isempty,xtmp) |...
                    cellfun(@isempty,ytmp) |...
                    cellfun(@isempty,ztmp));
    if ~isempty(ind_excl)
      fprintf('%s: WARNING: excluding entries with missing coordinates!\n',...
        mfilename);
      ind = setdiff([1:length(xtmp)],ind_excl);
      roi_info = roi_info(ind);
    end;

    % select points for this hemisphere
    if ismember(parms.hemilabel,fstrs)
      ind = find(strcmp({roi_info.(parms.hemilabel)},hemi));
    else
      % if no hemi column, use sign of x
      xvals = [roi_info.x];
      switch hemi
        case 'lh'
          ind = find(xvals <=  parms.hemi_thresh);
        case 'rh'
          ind = find(xvals >= -parms.hemi_thresh);
      end;
    end;
    roi_info = roi_info(ind);
    % create matrix of coordinates
    points = cat(2,[roi_info.x]',[roi_info.y]',[roi_info.z]');
    % create cell array of roinames
    if ismember(parms.namelabel,fstrs)
      roinames = regexprep({roi_info.(parms.namelabel)}',' ','_');
      roinames = cellfun(@(x) sprintf('ctx-%s-%s',hemi,x),...
                              roinames,'UniformOutput',false);
    else
      roinames = num2cell(1:length(roi_info))';
      roinames = cellfun(@(x) sprintf('ctx-%s-r%d',hemi,x),...
                              roinames,'UniformOutput',false);
    end;

    % get MNI flags for each point
    if ismember(parms.MNIlabel,fstrs)
      parms.MNIflag = [roi_info.(parms.MNIlabel)];
    end;

    mmil_mkdir(parms.outdir);

    % write roi_info for this hemisphere
    if parms.verbose
      fprintf('%s: writing roi_info to %s...\n',mfilename,fname_csv);
    end;
    mmil_struct2csv(roi_info,fname_csv);

    % create text file with roinames
    if parms.verbose
      fprintf('%s: writing roinames to %s...\n',...
        mfilename,fname_txt);
    end;
    fid = fopen(fname_txt,'wt');
    if fid==-1, error('failed to open %s for writing',fname_txt); end;
    fprintf(fid,'%s\n',roinames{:});
    fclose(fid);

    % sample points to surface
    if parms.verbose
      fprintf('%s: sampling points to surface...\n',mfilename);
    end;
    args = mmil_parms2args(parms,parms.surf_tags);
    vals = fs_tal2surf(points,'hemilist',hemi,args{:});

    % save weights to multi-frame mgh file
    if parms.verbose
      fprintf('%s: saving ROIs to %s...\n',mfilename,fname_out);
    end;
    vals = reshape(vals,[size(vals,1),1,1,size(vals,2)]);
    fs_save_mgh(vals,fname_out);
  end;
  fnames_weights{h} = fname_out;
  fnames_roinames{h} = fname_txt;
end;

