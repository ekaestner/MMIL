function rc_ready_retfit(subj,polstem,eccstem,varargin)
%function rc_ready_retfit(subj,polstem,eccstem,varargin)
%
% Purpose: copy retinotopy data and prepare for map fit with rc_retfit
%
% Required Parameters:
%   subj: FreeSurfer recon subject name
%   polstem: full file stem of polar angle retinotopy data (on cortical surface)
%     omit file ending (e.g. "_r-lh.mgh" or "_i-lh.mgh")
%     include full or relative path
%   eccstem: full file stem of eccentricity retinotopy data (on cortical surface)
%
% Optional Parameters:
%   'outdir': output directory
%     {default = 'retfit'}
%   'outstem': output file stem
%     {default = 'retfit'}
%   'roi_name': file stem of ROI file (e.g. lh.roi_name.label)
%     {default = 'v123'}
%   'r_min': minimum eccentricity (for stimulus presentation)
%     {default = 0.25}
%   'r_max': maximum eccentricity (for stimulus presentation)
%     {default = 15}
%   'logtrans_flag': [0|1] whether log transform was used for ecc stimulus
%     {default = 0}
%   'subjdir' - subjects directory (override SUBJECTS_DIR environment variable)
%     subjdir/subj should contain the FreeSurfer subject directory
%     {default = $SUBJECTS_DIR}
%   'map_v123_flag': [0|1] whether to model V1-V2-V3 complex or
%     a single area mapping entire hemifield
%     {default = 1}
%   'map_poly_flag': [0|1] use polynomial function to deform template
%     {default = 1}
%   'map_poly_order': order of polynomial function (n+1 additional parameters)
%      used to deform template
%     {default = 4}
%   'map_model_type': [0|1|2] model used for initial estimates of u and v
%     0: rectangle
%     1: wedge
%     2: radial wedge
%     {default = 2}
%   'map_area_name': area label if map_v123_flag=0
%     {default = 'v'}
%   'map_rev_polar_flag': [0|1] whether to reverse direction of polar angle
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output (including run_retfit.m)
%     {default = 0}
%
% Created:  02/14/11 by Don Hagler
% Last Mod: 12/05/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters
if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'outdir','retfit',[],...
  'outstem','retfit',[],...
  'roi_name','v123',[],...
  'r_min',0.25,[0,100],...
  'r_max',15,[1,100],...
  'logtrans_flag',false,[false true],...
  'subjdir',[],[],...
  'forceflag',false,[false true],...
... % retinotopy data / tksurfer
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'suffixlist',{'_r','_i'},{'_r','_i'},...
  'outstemlist',{'pol','ecc'},[],...
  'surf','sphere',{'white','pial','inflated','sphere'},...
  'smooth',10,[0,100],...
  'fthresh',0,[],...
  'fit_fthresh',0.01,[],...
  'fmid',1.5,[],...
  'fit_fmid',0.5,[],...
  'fslope',3,[],...
  'fit_fslope',3,[],...
  'revflag',false,[false true],...
  'sph_rot',{[45 0 90],[45 -20 -90]},[],...
... % retfit parameters
  'roi_dilate_niters',0,[0,1000],...
  'roi_rotation',65,[-180,180],...
  'roi_shift_u',0,[],...
  'roi_shift_v',0,[],...
  'roi_scale_u',0.7,[],...
  'roi_scale_v',0.7,[],...
  'prereg_nruns_quick',2,[],...
  'prereg_niter_quick',100,[],...
  'prereg_step_size_quick',[0.1,0.05,0.01],[],...
  'prereg_nruns',200,[],...
  'prereg_niter',100,[],...
  'map_v1_width_range',[1 1],[0.01,2],...
  'map_v2_width_range',[0.6 1],[0.01,2],...
  'map_v3_width_range',[0.5 1],[0.01,2],...
  'map_v1_length_range',[0.8 1.2],[0.1 10],...
  'map_v2_length_range',[0.8 1.2],[0.1 10],...
  'map_v3_length_range',[0.8 1.2],[0.1 10],...
  'map_poly_coef_range',[-5,5],[-100,100],...
  'map_wedge_fact_range',[0.7 1.3],[0.01,2],...
  'map_radial_wedge_fact_range',[0.05 0.4],[0,1],...
  'map_radial_offset_range',[1 4],[0,10],...
  'map_scale_u_range',[0.6 0.8],[0.01,1.0],...
  'map_scale_v_range',[0.3 0.6],[0.01,1.0],...
  'map_rotation_range',[-10 10],[-180,180],...
  'map_shift_u_range',[-0.1,0.1],[-1,1],...
  'map_shift_v_range',[-0.1,0.1],[-1,1],...
  'map_r_min_range',[1 1],[0,Inf],...
  'map_r_max_range',[12 12],[0,Inf],...
  'nruns',1,[],...
  'niter',2000,[],...
  'ecc_fact',1,[],...
  'smooth_fact',1,[],...
  'fold_fact',15,[],...
  'vacancy_fact',0,[0,Inf],...
  'max_outbound_penalty',20,[0 10000],...
  'data_smooth_sigma',0.1,[],...
  'map_v123_flag',true,[false true],...
  'map_poly_flag',true,[false true],...
  'map_poly_order',4,[1,10],...
  'map_model_type',2,[0,1,2],...
  'map_logtrans_flag',true,[false true],...
  'map_area_name','v',[],...
  'map_rev_polar_flag',false,[false true],...
  'cost_include_percentile',100,[],...
...
  'default_hemi','lh',[],...
});

parms.stemlist = {polstem,eccstem};
if mmil_isrelative(parms.outdir)
  parms.outdir = [pwd '/' parms.outdir];
end;

if isempty(parms.subjdir)
  parms.subjdir = deblank(getenv('SUBJECTS_DIR'));
  if isempty(parms.subjdir)
    error('Cannot find SUBJECTS_DIR environment variable');
  end
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% copy retinotopy data

data_outdir = [parms.outdir '/data'];
[succ,msg] = mkdir(data_outdir);
if ~succ, error('failed to create dir %s:\n%s',data_outdir,msg); end;

% check input files; if any missing, give error
for i=1:length(parms.stemlist)
  for j=1:length(parms.suffixlist)
    for h=1:length(parms.hemilist)
      fname_in = sprintf('%s%s-%s.mgh',parms.stemlist{i},...
        parms.suffixlist{j},parms.hemilist{h});
      if ~exist(fname_in,'file')
        error('file %s not found',fname_in);
      end;
    end;
  end;
end;

% check output files; if any missing, replace all
%   (this ensures log file accurately reflects source of data files)
copy_data_flag = parms.forceflag;
fname_log = [data_outdir '/data.log'];
if ~exist(fname_log,'file'), copy_data_flag = 1; end;
if ~copy_data_flag
  for i=1:length(parms.stemlist)
    for j=1:length(parms.suffixlist)
      for h=1:length(parms.hemilist)
        fname_out = sprintf('%s/%s%s-%s.mgh',data_outdir,...
          parms.outstemlist{i},parms.suffixlist{j},parms.hemilist{h});
        if ~exist(fname_out,'file'), copy_data_flag = 1; end;
      end;
    end;
  end;
end;

if copy_data_flag
  fid = fopen(fname_log,'wt');
  if fid<0
    error('failed to open %s for writing',fname_log);
  end;
  for i=1:length(parms.stemlist)
    for j=1:length(parms.suffixlist)
      for h=1:length(parms.hemilist)
        fname_in = sprintf('%s%s-%s.mgh',parms.stemlist{i},...
          parms.suffixlist{j},parms.hemilist{h});
        fname_out = sprintf('%s/%s%s-%s.mgh',data_outdir,...
          parms.outstemlist{i},parms.suffixlist{j},parms.hemilist{h});
        cmd = sprintf('cp %s %s',fname_in,fname_out);
        fprintf(fid,'%s\n',cmd);
        [status,result] = unix(cmd);
        if status
          fprintf(fid,'ERROR: %s\n',result);
          fclose(fid);
          error('cmd %s failed:\n%s',cmd,result);
        end;
      end;
    end;
  end;
  fclose(fid);
end;

parms.stemlist = parms.outstemlist;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create csh scripts for viewing retinoptopy data and fit

for i=1:length(parms.stemlist)
  fstem = parms.stemlist{i};
  indir = data_outdir;
  fname_view = sprintf('%s/view_%s.csh',parms.outdir,fstem);
  write_view_script(parms,fname_view,subj,fstem,indir)

  fstem = ['retfit_' parms.stemlist{i}];
  indir = [parms.outdir '/fit'];
  fname_view = sprintf('%s/view_%s.csh',parms.outdir,fstem);
  tmp_parms = parms;
  tmp_parms.fthresh = parms.fit_fthresh;
  tmp_parms.fmid = parms.fit_fmid;
  tmp_parms.fslope = parms.fit_fslope;
  write_view_script(tmp_parms,fname_view,subj,fstem,indir)
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create script for running retfit

fname_retfit = sprintf('%s/run_retfit.m',parms.outdir);
if ~exist(fname_retfit,'file') || parms.forceflag
  fid = fopen(fname_retfit,'wt');
  if fid==-1
    error('failed to open %s for writing',fname_retfit);
  end;
  fprintf(fid,'function run_retfit(hemi,test_run_flag,quick_flag,forceflag)\n');
  fprintf(fid,'%%function run_retfit([hemi],[test_run_flag],[quick_flag],[forceflag])\n');
  fprintf(fid,'%%\n');
  fprintf(fid,'%% Automatically generated script for running retinotopy fit\n');
  fprintf(fid,'%%\n');
  fprintf(fid,'%% Optional Parameters:\n');
  fprintf(fid,'%%   hemi: cortical hemisphere\n');
  fprintf(fid,'%%     {default = ''%s''}\n',parms.default_hemi);
  fprintf(fid,'%%   test_run_flag: [0|1] whether to make plots but not create final results\n');
  fprintf(fid,'%%     {default = 0}\n');
  fprintf(fid,'%%   quick_flag: [0|1] whether to do preregistration with multi-scale random search\n');
  fprintf(fid,'%%     {default = 0}\n');
  fprintf(fid,'%%   forceflag: [0|1] whether to overwrite existing output\n');
  fprintf(fid,'%%     {default = %d}\n',parms.forceflag);
  fprintf(fid,'%%\n');
  fprintf(fid,'%% Created:  %s\n',datestr(now,'mm/dd/yy'));
  fprintf(fid,'%% Last Mod: %s\n',datestr(now,'mm/dd/yy'));
  fprintf(fid,'%%\n');
  fprintf(fid,'\n');
  fprintf(fid,'if ~exist(''hemi'',''var'') | isempty(hemi), hemi = ''%s''; end;\n',parms.default_hemi);
  fprintf(fid,'if ~exist(''test_run_flag'',''var'') | isempty(test_run_flag), test_run_flag = 0; end;\n');
  fprintf(fid,'if ~exist(''quick_flag'',''var'') | isempty(quick_flag), quick_flag = 0; end;\n');
  fprintf(fid,'if ~exist(''forceflag'',''var'') | isempty(forceflag), forceflag = %d; end;\n',parms.forceflag);
  fprintf(fid,'\n');
  fprintf(fid,'subj = ''%s'';\n',subj);
  fprintf(fid,'\n');
  fprintf(fid,'parms = [];\n');
  fprintf(fid,'parms.subjdir = ''%s'';\n',parms.subjdir);
  fprintf(fid,'parms.polstem = ''%s/%s'';\n',data_outdir,parms.stemlist{1});
  fprintf(fid,'parms.eccstem = ''%s/%s'';\n',data_outdir,parms.stemlist{2});
  fprintf(fid,'parms.data_r_min = %0.4f;\n',parms.r_min);
  fprintf(fid,'parms.data_r_max = %0.4f;\n',parms.r_max);
  fprintf(fid,'parms.data_logtrans_flag = %d;\n',parms.logtrans_flag);
  fprintf(fid,'parms.roi_name = ''%s'';\n',parms.roi_name);
  fprintf(fid,'parms.outdir = ''%s'';\n',parms.outdir);
  fprintf(fid,'parms.outstem = ''%s'';\n',parms.outstem);
  fprintf(fid,'parms.hemilist = hemi;\n');
  fprintf(fid,'\n');
  fprintf(fid,'if ~test_run_flag\n');
  fprintf(fid,'  parms.save_results_flag = 1;\n');
  fprintf(fid,'  parms.plot_flag = 0;\n');
  fprintf(fid,'else\n');
  fprintf(fid,'  parms.save_results_flag = 0;\n');
  fprintf(fid,'  parms.plot_flag = 2;\n');
  fprintf(fid,'end;\n');
  fprintf(fid,'parms.forceflag = forceflag;\n');
  fprintf(fid,'\n');
  fprintf(fid,'%% ROI adjustment\n');
  fprintf(fid,'parms.roi_dilate_niters = %d;\n',parms.roi_dilate_niters);
  fprintf(fid,'switch hemi\n');
  fprintf(fid,'  case ''lh''\n');
  fprintf(fid,'    parms.roi_rotation = -%d;\n',parms.roi_rotation);
  fprintf(fid,'  case ''rh''\n');
  fprintf(fid,'    parms.roi_rotation = %d;\n',parms.roi_rotation);
  fprintf(fid,'end;\n');
  fprintf(fid,'parms.roi_shift_u = %0.4f;\n',parms.roi_shift_u);
  fprintf(fid,'parms.roi_shift_v = %0.4f;\n',parms.roi_shift_v);
  fprintf(fid,'parms.roi_scale_u = %0.4f;\n',parms.roi_scale_u);
  fprintf(fid,'parms.roi_scale_v = %0.4f;\n',parms.roi_scale_v);
  fprintf(fid,'\n');
  fprintf(fid,'%% pre-registration\n');
  fprintf(fid,'if quick_flag\n');
  fprintf(fid,'  parms.prereg_search_type = ''rand'';\n');
  fprintf(fid,'  parms.prereg_nruns = %d;\n',parms.prereg_nruns_quick);
  fprintf(fid,'  parms.prereg_niter = %d;\n',parms.prereg_niter_quick);
  fprintf(fid,'  parms.prereg_step_size = [%s];\n',...
    sprintf('%0.4f ',parms.prereg_step_size_quick));
  fprintf(fid,'else\n');
  fprintf(fid,'  parms.prereg_search_type = ''fmincon'';\n');
  fprintf(fid,'  parms.prereg_nruns = %d;\n',parms.prereg_nruns);
  fprintf(fid,'  parms.prereg_niter = %d;\n',parms.prereg_niter);
  fprintf(fid,'  parms.prereg_step_size = []; %% irrelevant\n');
  fprintf(fid,'end;\n');
  fprintf(fid,'parms.map_v123_flag = %d;\n',parms.map_v123_flag);
  fprintf(fid,'parms.map_poly_flag = %d;\n',parms.map_poly_flag);
  fprintf(fid,'parms.map_poly_order = %d;\n',parms.map_poly_order);
  fprintf(fid,'parms.map_model_type = %d;\n',parms.map_model_type);
  fprintf(fid,'parms.map_logtrans_flag = %d;\n',parms.map_logtrans_flag);
  fprintf(fid,'parms.map_area_name = ''%s'';\n',parms.map_area_name);
  fprintf(fid,'parms.map_rev_polar_flag = %d;\n',parms.map_rev_polar_flag);
  fprintf(fid,'parms.map_poly_coef_range = [%s];\n',...
    sprintf('%0.4f ',parms.map_poly_coef_range));
  fprintf(fid,'parms.map_radial_wedge_fact_range = [%s];\n',...
    sprintf('%0.4f ',parms.map_radial_wedge_fact_range));
  fprintf(fid,'parms.map_radial_offset_range = [%s];\n',...
    sprintf('%0.4f ',parms.map_radial_offset_range));
  fprintf(fid,'parms.map_scale_u_range = [%s];\n',...
    sprintf('%0.4f ',parms.map_scale_u_range));
  fprintf(fid,'parms.map_scale_v_range = [%s];\n',...
    sprintf('%0.4f ',parms.map_scale_v_range));
  fprintf(fid,'parms.map_rotation_range = [%s];\n',...
    sprintf('%0.4f ',parms.map_rotation_range));
  fprintf(fid,'parms.map_shift_u_range = [%s];\n',...
    sprintf('%0.4f ',parms.map_shift_u_range));
  fprintf(fid,'parms.map_shift_v_range = [%s];\n',...
    sprintf('%0.4f ',parms.map_shift_v_range));
  fprintf(fid,'parms.map_wedge_fact_range = [%s];\n',...
    sprintf('%0.4f ',parms.map_wedge_fact_range));
  fprintf(fid,'parms.map_r_min_range = [%s];\n',...
    sprintf('%0.4f ',parms.map_r_min_range));
  fprintf(fid,'parms.map_r_max_range = [%s];\n',...
    sprintf('%0.4f ',parms.map_r_max_range));
  fprintf(fid,'parms.map_v1_width_range = [%s];\n',...
    sprintf('%0.4f ',parms.map_v1_width_range));
  fprintf(fid,'parms.map_v2_width_range = [%s];\n',...
    sprintf('%0.4f ',parms.map_v2_width_range));
  fprintf(fid,'parms.map_v3_width_range = [%s];\n',...
    sprintf('%0.4f ',parms.map_v3_width_range));
  fprintf(fid,'parms.map_v1_length_range = [%s];\n',...
    sprintf('%0.4f ',parms.map_v1_length_range));
  fprintf(fid,'parms.map_v2_length_range = [%s];\n',...
    sprintf('%0.4f ',parms.map_v2_length_range));
  fprintf(fid,'parms.map_v3_length_range = [%s];\n',...
    sprintf('%0.4f ',parms.map_v3_length_range));
  fprintf(fid,'parms.ecc_fact = %0.4f;\n',parms.ecc_fact);
  fprintf(fid,'parms.vacancy_fact = %0.4f;\n',parms.vacancy_fact);
  fprintf(fid,'parms.max_outbound_penalty = %0.4f;\n',parms.max_outbound_penalty);
  fprintf(fid,'parms.data_smooth_sigma = %0.4f;\n',parms.data_smooth_sigma);
  fprintf(fid,'\n');
  fprintf(fid,'%% fine reg\n');
  fprintf(fid,'parms.nruns = %d;\n',parms.nruns);
  fprintf(fid,'parms.niter = %d;\n',parms.niter);
  fprintf(fid,'parms.smooth_fact = %0.4f;\n',parms.smooth_fact);
  fprintf(fid,'parms.fold_fact = %0.4f;\n',parms.fold_fact);
  fprintf(fid,'parms.cost_include_percentile = %0.4f;\n',parms.cost_include_percentile);
  fprintf(fid,'\n');
  fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n');
  fprintf(fid,'\n');
  fprintf(fid,'args = mmil_parms2args(parms);\n');
  fprintf(fid,'rc_retfit(subj,args{:})\n');
  fprintf(fid,'\n');
  fclose(fid);

  % change permissions so group can read and write
  cmd = sprintf('chmod ug+rw %s',fname_retfit);
  [status,result] = unix(cmd);
  if status
    fprintf('%s: WARNING: failed to set permissions for %s:\n%s\n',...
      mfilename,fname_retfit,result);
  end;
end;

return;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_view_script(parms,fname_out,subj,fstem,indir)
  tags = {'fname_out','fstem','indir','roi_name','label_dir','subjdir',...
    'forceflag','surf','smooth','fthresh','fmid','fslope','revflag',...
    'sph_rot','default_hemi'};

  parms.fstem = fstem;
  parms.fname_out = fname_out;
  parms.indir = indir;
  parms.label_dir = parms.outdir;
  
  args = mmil_parms2args(parms,tags);
  rc_write_viewscript(subj,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

