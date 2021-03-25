function [fname_wmparc,fname_clut] = mmil_fparc2wmparc(subj,varargin)
%function [fname_wmparc,fname_clut] = mmil_fparc2wmparc(subj,[options])
%
% Required parameters:
%   subj: FreeSurfer subject name
%
% Optional parameters:
%  'annotname': name of annotation
%    {default = 'fparc'}
%  'wmparcname' : output file stem
%    {default = 'wmparc'}
%  'outdir': output directory
%     if empty, will put in subjdir/subj/mri
%    {default = []}
%  'fnames_fparc': cell array of annotation files in fsaverage space
%    will be resampled to individual subject space before use
%    if supplied, fnames_aparc will be ignored
%    {default = []}
%  'fnames_aparc': cell array of annotation files (one for each hemisphere)
%    if empty, will use ?h.aparc.annot files in fspath/label
%    if not full path, assumed to be relative to fspath/label
%    {default = []}
%  'subjdir': directory containing FreeSurfer subject
%    {default = [getenv('FREESURFER_HOME') '/subjects']}
%  'verbose': [0|1] display status messages
%    {default = 1}
%  'forceflag': [0|1] overwrite existing output files
%    {default = 0}
%
% Output:
%   fname_wmparc: full path of output wmparc file
%   fname_clut: full path of output color lookup table text file
%
% Created:  09/10/15 by Don Hagler
% Last Mod: 09/13/15 by Don Hagler
%

fname_wmparc = [];
fname_clut = [];

if ~mmil_check_nargs(nargin,1), return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input parameters
parms = check_input(subj,varargin);

% set up output directory
parms = setup_outdir(parms);

% resample fparc from atlas to subject
if isempty(parms.fnames_aparc)
  parms = resample_fparc(parms);
end;

% create fparc+aseg from fparc
parms = create_fparc_seg(parms);

% create wmparc from fparc+aseg
parms = create_wmparc(parms);

% create color look-up table
parms = create_colorlut(parms);

% remove temporary directory
if parms.cleanup_flag
  cleanup_tmpdir(parms);
end;

fname_wmparc = parms.fname_wmparc;
fname_clut = parms.fname_clut;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input parameters
function parms = check_input(subj,options)
  parms = mmil_args2parms(options,{...
    'subj',subj,[],...
  ...
    'annotname','fparc',[],...
    'wmparcname','wmparc',[],...
    'outdir',[],[],...
    'fnames_fparc',[],[],...
    'fnames_aparc',[],[],...
    'subjdir',[getenv('FREESURFER_HOME') '/subjects'],[],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ...
    'tmpdir','tmp_fparc2wmparc',[],...
    'cleanup_flag',false,[false true],...
    'verbose',false,[false true],...
    'hemilist',{'lh','rh'},{'lh','rh'},...
    'gm_codes',[1000,2000],[],...
    'wm_codes',[3000,4000],[],...
    'extra_lut_lines',{...
      '0   Unknown                                 0   0   0   0'...
      '1   Left-Cerebral-Exterior                  70  130 180 0'...
      '2   Left-Cerebral-White-Matter              245 245 245 0'...
      '3   Left-Cerebral-Cortex                    205 62  78  0'...
      '4   Left-Lateral-Ventricle                  120 18  134 0'...
      '5   Left-Inf-Lat-Vent                       196 58  250 0'...
      '6   Left-Cerebellum-Exterior                0   148 0   0'...
      '7   Left-Cerebellum-White-Matter            220 248 164 0'...
      '8   Left-Cerebellum-Cortex                  230 148 34  0'...
      '9   Left-Thalamus                           0   118 14  0'...
      '10  Left-Thalamus-Proper                    0   118 14  0'...
      '11  Left-Caudate                            122 186 220 0'...
      '12  Left-Putamen                            236 13  176 0'...
      '13  Left-Pallidum                           12  48  255 0'...
      '14  3rd-Ventricle                           204 182 142 0'...
      '15  4th-Ventricle                           42  204 164 0'...
      '16  Brain-Stem                              119 159 176 0'...
      '17  Left-Hippocampus                        220 216 20  0'...
      '18  Left-Amygdala                           103 255 255 0'...
      '19  Left-Insula                             80  196 98  0'...
      '20  Left-Operculum                          60  58  210 0'...
      '21  Line-1                                  60  58  210 0'...
      '22  Line-2                                  60  58  210 0'...
      '23  Line-3                                  60  58  210 0'...
      '24  CSF                                     60  60  60  0'...
      '25  Left-Lesion                             255 165 0   0'...
      '26  Left-Accumbens-area                     255 165 0   0'...
      '27  Left-Substancia-Nigra                   0   255 127 0'...
      '28  Left-VentralDC                          165 42  42  0'...
      '29  Left-undetermined                       135 206 235 0'...
      '30  Left-vessel                             160 32  240 0'...
      '31  Left-choroid-plexus                     0   200 200 0'...
      '32  Left-F3orb                              100 50  100 0'...
      '33  Left-lOg                                135 50  74  0'...
      '34  Left-aOg                                122 135 50  0'...
      '35  Left-mOg                                51  50  135 0'...
      '36  Left-pOg                                74  155 60  0'...
      '37  Left-Stellate                           120 62  43  0'...
      '38  Left-Porg                               74  155 60  0'...
      '39  Left-Aorg                               122 135 50  0'...
      '40  Right-Cerebral-Exterior                 70  130 180 0'...
      '41  Right-Cerebral-White-Matter             0   225 0   0'...
      '42  Right-Cerebral-Cortex                   205 62  78  0'...
      '43  Right-Lateral-Ventricle                 120 18  134 0'...
      '44  Right-Inf-Lat-Vent                      196 58  250 0'...
      '45  Right-Cerebellum-Exterior               0   148 0   0'...
      '46  Right-Cerebellum-White-Matter           220 248 164 0'...
      '47  Right-Cerebellum-Cortex                 230 148 34  0'...
      '48  Right-Thalamus                          0   118 14  0'...
      '49  Right-Thalamus-Proper                   0   118 14  0'...
      '50  Right-Caudate                           122 186 220 0'...
      '51  Right-Putamen                           236 13  176 0'...
      '52  Right-Pallidum                          13  48  255 0'...
      '53  Right-Hippocampus                       220 216 20  0'...
      '54  Right-Amygdala                          103 255 255 0'...
      '55  Right-Insula                            80  196 98  0'...
      '56  Right-Operculum                         60  58  210 0'...
      '57  Right-Lesion                            255 165 0   0'...
      '58  Right-Accumbens-area                    255 165 0   0'...
      '59  Right-Substancia-Nigra                  0   255 127 0'...
      '60  Right-VentralDC                         165 42  42  0'...
      '61  Right-undetermined                      135 206 235 0'...
      '62  Right-vessel                            160 32  240 0'...
      '63  Right-choroid-plexus                    0   200 221 0'...
      '64  Right-F3orb                             100 50  100 0'...
      '65  Right-lOg                               135 50  74  0'...
      '66  Right-aOg                               122 135 50  0'...
      '67  Right-mOg                               51  50  135 0'...
      '68  Right-pOg                               74  155 60  0'...
      '69  Right-Stellate                          120 62  43  0'...
      '70  Right-Porg                              74  155 60  0'...
      '71  Right-Aorg                              122 135 50  0'...
      '72  5th-Ventricle                           120 190 150 0'...
      '73  Left-Interior                           122 135 50  0'...
      '74  Right-Interior                          122 135 50  0'...
      '77  WM-hypointensities                      200 70  255 0'...
      '78  Left-WM-hypointensities                 255 148 10  0'...
      '79  Right-WM-hypointensities                255 148 10  0'...
      '80  non-WM-hypointensities                  164 108 226 0'...
      '81  Left-non-WM-hypointensities             164 108 226 0'...
      '82  Right-non-WM-hypointensities            164 108 226 0'...
      '83  Left-F1                                 255 218 185 0'...
      '84  Right-F1                                255 218 185 0'...
      '85  Optic-Chiasm                            234 169 30  0'...
      '192 Corpus_Callosum                         250 255 50  0'...
      '5001 Left-UnsegmentedWhiteMatter            20  30  40  0'...
      '5002 Right-UnsegmentedWhiteMatter           20  30  40  0'...
      },[],...
  });

  parms.nhemi = length(parms.hemilist);

  if isempty(parms.subjdir)
    parms.subjdir = getenv('SUBJECTS_DIR');
    if isempty(parms.subjdir)
      error('SUBJECTS_DIR not defined as an environment variable');
    end;
  else
    setenv('SUBJECTS_DIR',parms.subjdir);
  end;

  % set fspath
  parms.fspath = sprintf('%s/%s',parms.subjdir,parms.subj);
  if ~exist(parms.fspath,'file')
    error('FS recon dir %s not found',parms.fspath);
  end;

  % check fparc file
  if ~isempty(parms.fnames_fparc)
    if ~iscell(parms.fnames_fparc)
      parms.fnames_fparc = {parms.fnames_fparc};
    end;
    for f=1:length(parms.fnames_fparc)
      if ~exist(parms.fnames_fparc{f},'file')
        error('fparc annot file %s not found',parms.fnames_fparc{f});
      end;
    end;
  end;

  % check aparc file if no fparc file
  if isempty(parms.fnames_fparc)
    if isempty(parms.fnames_aparc)
      for h=1:parms.nhemi
        hemi = parms.hemilist{h};
        parms.fnames_aparc{h} = sprintf('%s/label/%s.aparc.annot',...
          parms.fspath,hemi);
      end;
    else
      if ~iscell(parms.fnames_aparc)
        parms.fnames_aparc = {parms.fnames_aparc};
      end;
      if length(parms.fnames_aparc) ~= parms.nhemi
        error('must have %d elements in fnames_aparc (have %d)',...
          parms.nhemi,length(parms.fnames_aparc));
      end;
    end;
    for h=1:parms.nhemi
      hemi = parms.hemilist{h};
      if mmil_isrelative(parms.fnames_aparc{h})
        parms.fnames_aparc{h} = [parms.fspath '/label/' parms.fnames_aparc{h}];
      end;
      if ~exist(parms.fnames_aparc{h},'file')
        error('file %s not found',parms.fnames_aparc{h});
      end;
      parms.aparc_hemis{h} = hemi;
      parms.aparc_names{h} = 'aparc';
    end;
  else
    parms.fnames_aparc = [];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create temporary output directory
function parms = setup_tmpdir(parms)
  % set tmpdir to be relative to outdir
  if mmil_isrelative(parms.tmpdir)
    parms.tmpdir = [parms.outdir '/' parms.tmpdir];
  end;
  mmil_mkdir(parms.tmpdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up output directory
function parms = setup_outdir(parms)
  if isempty(parms.outdir)
    parms.outdir = [parms.fspath '/mri'];
    parms = setup_tmpdir(parms);
  else
    parms = setup_tmpdir(parms);
    % create copy of subj dir with links for most contents
    %   but writable copies for mri and label subdirs
    subjdir_new = parms.tmpdir;
    fspath_new = [subjdir_new '/' parms.subj];
    % remove existing dir if already exists and forceflag
    if exist(fspath_new,'dir') && parms.forceflag
      cmd = sprintf('rm -r %s',fspath_new);
      [s,r] = unix(cmd);
      if s, error('cmd failed: %s:\n%s',cmd,r); end;
    end;
    if ~exist(fspath_new,'dir')
      if parms.verbose
        fprintf('%s: creating writable copy of %s...\n',mfilename,parms.subj);
      end;
      % create new subj dir
      mmil_mkdir(fspath_new);
      % create links to all content
      cmd = sprintf('cd %s; ln -s %s/* .',fspath_new,parms.fspath);
      [s,r] = unix(cmd); if s, error('cmd failed: %s:\n%s',cmd,r); end;
      % remove links for mri and label
      cmd = sprintf('rm %s/label %s/mri',fspath_new,fspath_new);
      [s,r] = unix(cmd); if s, error('cmd failed: %s:\n%s',cmd,r); end;
      % create writable copies of mri and label subdirs
      cmd = sprintf('cp -r %s/label %s/mri %s',...
        parms.fspath,parms.fspath,fspath_new);
      [s,r] = unix(cmd); if s, error('cmd failed: %s:\n%s',cmd,r); end;
    end;
    % reset subjdir, subj, fspath
    parms.subjdir = subjdir_new;
    parms.fspath = fspath_new;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% resample fparc annotation file from fsaverage to individual
function parms = resample_fparc(parms)
  for f=1:length(parms.fnames_fparc)
    fname_in = parms.fnames_fparc{f};
    if parms.verbose
      fprintf('%s: resampling annotation file %s from fsaverage to %s...\n',...
        mfilename,fname_in,parms.subj);
    end;
    % call fs_annot2annot
    tmp_parms = [];
    tmp_parms.outdir = [parms.fspath '/label'];
    tmp_parms.source_subj = 'fsaverage';
    tmp_parms.subj = parms.subj;
    tmp_parms.subjdir = parms.subjdir;
    tmp_parms.verbose = parms.verbose;
    tmp_parms.forceflag = parms.forceflag;
    args = mmil_parms2args(tmp_parms);
    fname_out = fs_annot2annot(fname_in,args{:});
    n = regexp(fname_out,'(?<hemi>[lr]h)\.(?<name>.+)\.annot$','names');
    if isempty(n)
      error('unexpected naming convention for aparc file %s\n',fname_out);
    end;
    parms.fnames_aparc{f} = fname_out;
    parms.aparc_hemis{f} = n.hemi;
    parms.aparc_names{f} = n.name;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create fparc+aseg from fparc
function parms = create_fparc_seg(parms)
  if parms.verbose
    fprintf('%s: creating fparc+aseg...\n',mfilename);
  end;
  fname_out = sprintf('%s/%s+aseg.mgz',parms.outdir,parms.annotname);
  if ~exist(fname_out,'file') || parms.forceflag
    fname_tmp = sprintf('%s/mri/%s+aseg.mgz',parms.fspath,parms.annotname);
    if ~exist(fname_tmp,'file') || parms.forceflag
      cmd = [];
      cmd = sprintf('%s setenv SUBJECTS_DIR %s\n',cmd,parms.subjdir);
      cmd = sprintf('%s mri_aparc2aseg --s %s',cmd,parms.subj);
      cmd = sprintf('%s   --annot %s',cmd,parms.annotname);
      cmd = sprintf('%s   --annot-table %s/label/%s.ctab',cmd,parms.fspath,parms.annotname);
      cmd = sprintf('%s   --volmask',cmd);
      cmd = sprintf('%s   --o %s\n',cmd,fname_tmp);
      [s,r] = unix(cmd);
      if s || ~exist(fname_tmp,'file')
        error('cmd failed: %s:\n%s',cmd,r);
      end;
    end;
    if ~strcmp(fname_tmp,fname_out)
      cmd = sprintf('cp -p %s %s',fname_tmp,fname_out);
      [s,r] = unix(cmd);
      if s || ~exist(fname_out,'file')
        error('cmd failed: %s:\n%s',cmd,r);
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create wmparc from fparc+aseg
function parms = create_wmparc(parms)
  if parms.verbose
    fprintf('%s: creating wmparc...\n',mfilename);
  end;
  fname_out = sprintf('%s/%s.mgz',parms.outdir,parms.wmparcname);
  if ~exist(fname_out,'file') || parms.forceflag
    cmd = [];
    cmd = sprintf('%s setenv SUBJECTS_DIR %s\n',cmd,parms.subjdir);
    cmd = sprintf('%s mri_aparc2aseg --s %s',cmd,parms.subj);
    cmd = sprintf('%s   --annot %s',cmd,parms.annotname);
    cmd = sprintf('%s   --annot-table %s/label/%s.ctab',cmd,parms.fspath,parms.annotname);
    cmd = sprintf('%s   --wmparc-dmax 5',cmd);
    cmd = sprintf('%s   --labelwm --hypo-as-wm --rip-unknown',cmd);
    cmd = sprintf('%s   --volmask',cmd);
    cmd = sprintf('%s   --ctxseg %s+aseg.mgz',cmd,parms.annotname);
    cmd = sprintf('%s   --o %s',cmd,fname_out);
    [s,r] = unix(cmd);
    if s || ~exist(fname_out,'file')
      error('cmd failed: %s:\n%s',cmd,r);
    end;
  end;
  parms.fname_wmparc = fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create color look-up table
function parms = create_colorlut(parms)
  parms.fname_ctab = sprintf('%s/label/%s.ctab',...
                             parms.fspath,parms.annotname);
  parms.fname_clut = sprintf('%s/%s_colorLUT.txt',...
                             parms.outdir,parms.wmparcname);
  if ~exist(parms.fname_clut,'file') || parms.forceflag
    % load ctab file
    [roicodes, roinames, rgbv] = fs_colorlut(parms.fname_ctab);
    nroi = length(roicodes);

    % set new roicodes, roinames, and rgbv for wmparc ROIs
    nnew = 2*parms.nhemi*nroi;
    new_roicodes = zeros(nnew,1);
    new_roinames = cell(nnew,1);
    new_rgbv = zeros(nnew,4);
    k = 1;
    for h=1:parms.nhemi
      hemi = parms.hemilist{h};
      base_code = parms.gm_codes(h);
      for r=1:nroi
        new_roicodes(k) = base_code + roicodes(r);
        new_roinames{k} = sprintf('ctx-%s-%s',hemi,roinames{r});
        new_rgbv(k,:) = rgbv(r,:);
        k = k + 1;
      end;
    end;
    for h=1:parms.nhemi
      hemi = parms.hemilist{h};
      base_code = parms.wm_codes(h); 
      for r=1:nroi
        new_roicodes(k) = base_code + roicodes(r);
        new_roinames{k} = sprintf('wm-%s-%s',hemi,roinames{r});
        new_rgbv(k,:) = rgbv(r,:);
        k = k + 1;
      end;
    end;

    % create color LUT file
    fid = fopen(parms.fname_clut,'wt');
    if fid<0, error('failed to open %s for writing',parms.fname_clut); end;
    % write extra lines
    for i=1:length(parms.extra_lut_lines)
      fprintf(fid,'%s\n',parms.extra_lut_lines{i});
    end;
    % write codes, names, and colors for new ROIs
    for r=1:nnew
      rgb = sprintf(' %3d',new_rgbv(r,:));
      fprintf(fid,'%4d  %-40s %s\n',new_roicodes(r),new_roinames{r},rgb);
    end;
    fclose(fid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove temporary directory
function cleanup_tmpdir(parms)
  cmd = sprintf('rm -r %s',parms.tmpdir);
  [s,r] = unix(cmd);
  if s, error('cmd failed: %s:\n%s',cmd,r); end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

