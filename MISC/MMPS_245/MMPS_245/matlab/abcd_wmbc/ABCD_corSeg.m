function [segStruct regStruct] = corSeg_amd(varargin)
% Perform segmentation on a T1-weighted volume (with COR or DICOM inputs)
%
% segStruct = corSeg(varargin);
%
% Inputs: (in any order)
% -- imgDir: a directory with T1w COR or DICOM files   (DEFAULT: 'orig')
% -- 'reslice on/off': relice to 1mm 256^3 coronal view in RB atlas space
%    (DEFAULT: off)
% -- 'headn3 on/off': N3 bias-field correction for head     (DEFAULT: on)
% -- 'headn3more on/off': N3 bias-field correction for head 4 times more   (DEFAULT: false)
% -- 'brainn3 on/off': N3 bias-field correction for brain   (DEFAULT: on)
% -- 'brainn3more on/off': N3 bias-field correction for brain 4 times more   (DEFAULT: false)
% -- 'histeqatl on/off': histogram eqalization atlas to vol (DEFAULT: false)
% -- 'histeqatlbf on/off': histogram eqalization atlas to vol +bias field  (DEFAULT: false)
% -- 'histeqatlmuk on/off': histogram eqalization atlas to vol by using Muk's  (DEFAULT: false)
% -- 'histeqvol on/off': histogram eqalization vol to atlas (DEFAULT: false)
% -- 'histeqvolmuk on/off': histogram eqalization vol to atlas by using Muk's (DEFAULT: false)
% -- 'mesh on/off': LevelSet to clean up seg           (DEFAULT: on)
% -- 'save on/off': save segStruct.mat in the cwd      (DEFAULT: on)
% -- strID: a string starting with '_' to ID the saves (DEFAULT: '')
% -- 'display on/off': call showSeg at the end         (DEFAULT: off)
% -- feature_atlas:                (DEFAULT: feature_buckner_mpr_256)

% ------------------------------------------------------------ 
% Don't forget: Need to check orig1's Mvxl2lph
%
% orig1.Mvxl2lph = origa.Mvxl2lph; % !@#
% orig1 = vol_add_fields(orig1);
% ------------------------------------------------------------ 

fp = 1;
totTime = 0;

DO_SAVE = false;
DISPLAY = false;
DO_HEAD_N3 = true;
DO_HEAD_N3_MORE = false;
DO_BRAIN_N3 = true;
DO_BRAIN_N3_MORE = false;
DO_HISTEQ_ATL = false;
DO_HISTEQ_ATLBF = false;
DO_HISTEQ_ATLMUK = false;
DO_HISTEQ_VOL = false;
DO_HISTEQ_VOLMUK = false;
USE_MGHBRAIN = false;
DO_RESLICE = false;
DO_MESH = true;
DO_WRITEDICOM = false;
DO_WRITESCDCM = false;
DO_WRITEJPEG = false;
%USE_1MM = false;
USE_1MM = true;
WM_label = 1;

% ------------------------------------------------------------ 
% Inputs
% ------------------------------------------------------------ 

imgDir = 'orig';
strID = '';
atlas_ext = '';

% Inputs in any order
for ii = 1 : nargin
  if ischar(varargin{ii})
    switch varargin{ii}
     case 'atlas edited',
      atlas_ext = '_edited';
     case '1mm on',
      USE_1MM = true;
     case '1mm off',
      USE_1MM = false;    
     case 'jpeg on',
      DO_WRITEJPEG = true;
     case 'jpeg off',
      DO_WRITEJPEG = false; 
     case 'scdcms on',
      DO_WRITESCDCM = true;
     case 'scdcms off',
      DO_WRITESCDCM = false; 
     case 'dcms on',
      DO_WRITEDICOM = true;
     case 'dcms off',
      DO_WRITEDICOM = false;
     case 'headn3 on',
      DO_HEAD_N3 = true;
     case 'headn3 off',
      DO_HEAD_N3 = false;
     case 'headn3more on',
      DO_HEAD_N3_MORE = true;
     case 'headn3more off',
      DO_HEAD_N3_MORE = false;      
     case 'brainn3 on',
      DO_BRAIN_N3 = true;
     case 'brainn3 off',
      DO_BRAIN_N3 = false;
     case 'brainn3more on',
      DO_BRAIN_N3_MORE = true;
     case 'brainn3more off',
      DO_BRAIN_N3_MORE = false;
     case 'histeqatl on',
      DO_HISTEQ_ATL = true;
     case 'histeqatl off',
      DO_HISTEQ_ATL = false; 
     case 'histeqatlmuk on',
      DO_HISTEQ_ATLMUK = true;
     case 'histeqatlmuk off',
      DO_HISTEQ_ATLMUK = false;    
     case 'histeqatlbf on',
      DO_HISTEQ_ATLBF = true;
     case 'histeqatlbf off',
      DO_HISTEQ_ATLBF = false; 
     case 'histeqvol on',
      DO_HISTEQ_VOL = true;
     case 'histeqvol off',
      DO_HISTEQ_VOL = false;
     case 'histeqvolmuk on',
      DO_HISTEQ_VOLMUK = true;
     case 'histeqvolmuk off',
      DO_HISTEQ_VOLMUK = false;   
     case 'mghbrain on',
      USE_MGHBRAIN = true;
     case 'mghbrain off',
      USE_MGHBRAIN = false;    
     case 'reslice on',
      DO_RESLICE = true;
     case 'reslice off',
      DO_RESLICE = false;  
     case 'mesh on',
      DO_MESH = true;
     case 'mesh off',
      DO_MESH = false;      
     case 'save on', 
      DO_SAVE = true;
     case 'save off', 
      DO_SAVE = false;
     case 'display on', 
      DISPLAY = true;
     case 'display off', 
      DISPLAY = false;
     otherwise,
      if strncmp(varargin{ii}, '_', 1)
        strID = varargin{ii};
      else
        imgDir = varargin{ii};
      end
    end
  elseif isstruct(varargin{ii})
    if isfield(varargin{ii}, 'mu_1')
      feature_atlas = varargin{ii};
    elseif isfield(varargin{ii}, 'mufact') % Allow adjustment of expected image intensities for each label (AMD)
      feature_atlas_adjustment = varargin{ii};
    elseif isfield(varargin{ii}, 'imgs') % Allow dircet input of image volume (AMD)
      imgDir = varargin{ii};
    end
  end
end

if (DO_HEAD_N3_MORE)
    DO_HEAD_N3=true;
end

if (DO_BRAIN_N3_MORE)
    DO_BRAIN_N3=true;
end
% Can only have one option for contrast adjustment
if (DO_HISTEQ_ATL)
    DO_HISTEQ_VOL = false;
    DO_HISTEQ_ATLBF = false;
    DO_HISTEQ_VOLMUK =false;
    DO_HISTEQ_ATLMUK =false;
end

if (DO_HISTEQ_VOL)
    DO_HISTEQ_ATL = false;
    DO_HISTEQ_ATLBF = false;
    DO_HISTEQ_VOLMUK =false;
    DO_HISTEQ_ATLMUK =false;
end

if (DO_HISTEQ_VOLMUK)
    DO_HISTEQ_ATL = false;
    DO_HISTEQ_ATLBF = false;
    DO_HISTEQ_VOL =false;
    DO_HISTEQ_ATLMUK =false;
end

if (DO_HISTEQ_ATLBF)
    DO_HISTEQ_VOL = false;
    DO_HISTEQ_ATL = false;
    DO_HISTEQ_VOLMUK =false;
    DO_HISTEQ_ATLMUK =false;
end

if (DO_HISTEQ_ATLMUK)
    DO_HISTEQ_VOL = false;
    DO_HISTEQ_ATL = false;
    DO_HISTEQ_ATLBF = false;
    DO_HISTEQ_VOLMUK =false;
end

% ------------------------------------------------------------  
% Read in T1w volume
% ------------------------------------------------------------ 

tic;  

if isstruct(imgDir)
  fprintf(fp, 'Using input volume ... ');
  orig = imgDir;
  orig.sf = 1.;
  DO_SAVE = false;
else

fprintf(fp, 'Read in %s ... ', imgDir);
if~(isdir(imgDir))
%    [pathstr, fname, fext, versn] = fileparts(imgDir);
    [pathstr, fname, fext] = fileparts(imgDir);
    if (strcmp(fext, '.mgh'))
        orig = read_mgh(imgDir);
        orig.sf = 1;
    elseif (strcmp(fext, '.mat'))
        t = load(imgDir);
        orig = getfield(t, fname);
        orig.sf=1;
        clear t;
    else
        fprintf('\nFile %s doesn''t exist or format is not correct\n\n', imgDir);
        error('Exiting ...');
    end
        
else
    if ~exist(imgDir, 'dir');
        fprintf('\nDirectory %s doesn''t exist\n\n', imgDir);
        error('Exiting ...');
    end

    cd(imgDir);
    if exist('COR-.info', 'file');
        COR_DIR = true;
    else
        COR_DIR = false;
    end
    cd('..');

    if COR_DIR

        orig = read_cor(imgDir);
        orig.imgs = orig.imgs;
        orig.sf = 1.;
        fprintf('COR files ... ');

    else % assume we have DICOMs
        dcmdir = [imgDir '//dcm//'];
        if (exist(dcmdir, 'dir'))
            bdcm = true;
        else
            bdcm = false;
        end
        try
            if (bdcm)
                orig = read_dicomdir(dcmdir, true);
            else
                orig = read_dicomdir(imgDir, true);
            end
            dcmhdr = orig.dcminfo;
            segStruct.dcmhdr = orig.dcminfo;
            cd(imgDir);
            fprintf('DICOM files ... ');
        catch
            fprintf('\nUnable to read DICOMs from directory %s\n\n', imgDir);
            error('Exiting ...');
        end
        
    end
end

end

% Rescaling if needed
if ~isfield(orig,'maxI'), [orig.minI orig.maxI] = deal(0,max(orig.imgs(:))); end % AMD

if orig.maxI > 4096
    orig.imgs = orig.imgs*4096/(orig.maxI-orig.minI)+orig.minI;
    orig.maxI=4096;
    orig.minI=0;
end

t = toc; 
totTime = totTime + t;
fprintf(fp, 'done (%.1f s)\n', t);
% ------------------------------------------------------------  
% Unbias the volumes (N3)   for head                                   
% ------------------------------------------------------------  

if DO_HEAD_N3
    
fprintf(fp, 'Unbias HEAD ... ');
tic;  

estop = 0.001;
maxiter = 100;

orig = vol_correct_bias_field_n3(orig, estop, maxiter);

if (DO_HEAD_N3_MORE)
    for iter=1:4
        orig = vol_correct_bias_field_n3(orig, estop, maxiter, 200, 0.15, 0.01, 50, 4);
    end
end


t = toc; 
totTime = totTime + t;
fprintf(fp, 'done (%.1f s)\n', t);

else
    fprintf(fp, 'Skip N3 for head !! \n');
end

% % ------------------------------------------------------------
% Do the Skull Stripping 
% ------------------------------------------------------------

fprintf(fp, 'Skull stripping ... ');
tic;
sampling =[4 4 4];
 nK = [5 5 5]; % Gen says (Apr 5, 2006) 5 is best
 boutput = false;
 sf = 1;
[bm_gen, M_hatl_to_head_rb, regStruct] = dctMorph_SkullStrip(orig, sampling, nK, boutput, sf, false);
segStruct.brainmesh = SkullStripAMDwrapper(orig, regStruct);
segStruct.MI = regStruct.min_cost_rb;
orig.imgs = orig.imgs *regStruct.sfm;
orig.maxI = 600;
orig.minI= 0;
clear regStruct volm;
%segStruct.brainmesh =  skullstrip2006(orig, bm_gen, M_hatl_to_head_rb);
orig_brain = getmaskvol(orig, segStruct.brainmesh, eye(4,4));
t = toc;
totTime = totTime + t;
fprintf(fp, 'done (%.1f s)\n', t);

% ------------------------------------------------------------  
% Unbias the volumes (N3) for brain                                  
% ------------------------------------------------------------  

if DO_BRAIN_N3
    
fprintf(fp, 'Unbias BRAIN ... ');
tic;  

estop = 0.001;
maxiter = 100;

orig_brain = vol_correct_bias_field_n3(orig_brain, estop, maxiter);

if (DO_BRAIN_N3_MORE)
    for iter=1:2
        orig_brain = vol_correct_bias_field_n3(orig_brain, estop, maxiter);
    end
end

t = toc;
totTime = totTime + t;
fprintf(fp, 'done (%.1f s)\n', t);

else
    fprintf(fp, 'Skip N3 for brain !! \n');
end

% ------------------------------------------------------------
% Do alignment(rb+affine+nm) to atlas
% ------------------------------------------------------------
fprintf(fp, 'Alignment to atlas and creat voxel space mapping between atlas and subject... ');
tic;

[regStruct, M_atl_to_vol_rb, volm] = brain_reg2atl(orig_brain, M_hatl_to_head_rb, DO_RESLICE);
segStruct.sfm = regStruct.sfm;
if (DO_RESLICE)
    orig =  vol_reslice(orig, 0, 0, [256 256 256 1 1 1], 2, M_atl_to_vol_rb);
    orig_brain =  vol_reslice(orig_brain, 0, 0, [256 256 256 1 1 1], 2, M_atl_to_vol_rb);
    lph = segStruct.brainmesh.vertices;
    lph(:,4)=1;
    lph_rbatl = (inv(M_atl_to_vol_rb)*lph')';
    segStruct.brainmesh.vertices = lph_rbatl(:,1:3);
    clear lph lph_rbatl;
    M_atl_to_vol_rb =eye(4,4);
end
vxlmap = getDCTVxlMapping(orig_brain, regStruct, volm);
t=toc;
totTime = totTime + t;
fprintf(fp, 'done (%.1f s)\n', t)

% ------------------------------------------------------------
% Save regStruct and vxlmap
% ------------------------------------------------------------

if DO_SAVE

    fprintf(fp, 'Saving %s/regStruct%s.mat ... ', pwd, strID);
    tic;

    save(['regStruct' strID '.mat'], 'regStruct');

    t = toc;
    totTime = totTime + t;
    fprintf(fp, 'done (%.1f s)\n', t);
    
    fprintf(fp, 'Saving %s/vxlmap%s.mat ... ', pwd, strID);
    tic;

%    save(['vxlmap' strID '.mat'], 'vxlmap');

    t = toc;
    totTime = totTime + t;
    fprintf(fp, 'done (%.1f s)\n', t);
else

    fprintf(fp, 'SKIPPING saving to regStruct.mat ... \n');

end

%clear regStruct volm;
clear volm;
% ------------------------------------------------------------
% Run estimateTissueStats_atl to correct fa AND globI
% ------------------------------------------------------------

fprintf(fp, 'Correcting contrast for atlas and volume... ');
tic;
orig_brain.imgs = orig_brain.imgs * segStruct.sfm;
orig_brain.maxI=600;
orig_brain.minI =0;
if (USE_1MM)
    load(sprintf('ctxatl_1_4_rbm%s.mat',atlas_ext));
else
    load(sprintf('ctxatl_4_4_rbm%s.mat',atlas_ext));
end
if (DO_HISTEQ_ATL)
    load hpvol_rbm;
    brain_atl = getVolFromVxlMap(orig_brain , vxlmap,  0, 2);
    ctxatl.locstats = adjustCTXfa(ctxatl.locstats, hpvol, brain_atl, 0.5, 10, 60);
    clear hpvol brain_atl;
end
% change volume's intensity directly
% if (DO_HISTEQ_VOL)   
%     load hpvol;
%     volm_subj = getVolFromVxlMap(hpvol , vxlmap,  1, 2);
%     orig_brain = adjustContrast_ctxvol(orig_brain, volm_subj);
%     clear hpvol volm_subj;
% end
% 
% if (DO_HISTEQ_ATLBF)
%     NUM_SMOOTHS = 10;
%     ALPHA = 1;
%     % Compute the bias field
%     load hpvol;
%     volm_subj = getVolFromVxlMap(hpvol , vxlmap,  1, 2);
%     clear hpvol;
%     brain = adjustContrast_ctxvol(orig_brain, volm_subj);
%     volSR = volSmoothRatio(volm_subj, brain, NUM_SMOOTHS, ALPHA);
%     orig_brain.imgs = orig_brain.imgs .* volSR.imgs;
%     clear brain volSR;
%     ctxatl.locstats = adjustContrast_ctxfa( ctxatl.locstats, orig_brain, volm_subj);
%     clear volm_subj;
% end


t = toc;
totTime = totTime + t;
fprintf(fp, 'done (%.1f s)\n', t);

fprintf(fp, 'Do the voxel Labeling by MRF ... ');
tic;

if exist('feature_atlas_adjustment','var') % (AMD)
  ctxatl.feature_atlas_adjustment = feature_atlas_adjustment;
end

% Temporary hack (AMD)
%save('snap.mat'); % Save snapshot
% load('snap.mat');

segStruct.topLabel = ctxSeg(orig_brain , ctxatl, vxlmap, 3, 5); % AMD: what's does the last argument do? 2 vs. 5 iterations?
segStruct.Left_cm3 = segStruct.topLabel.Left_cm3;
segStruct.Right_cm3 = segStruct.topLabel.Right_cm3;

clear ctxatl vxlmap;
t = toc;
totTime = totTime + t;
fprintf(fp, 'done (%.1f s)\n', t);

% ------------------------------------------------------------
% Save segStruct
% ------------------------------------------------------------

segStruct.topLabel.imgs = double(segStruct.topLabel.imgs);
segStruct.vol2 = orig_brain;
clear orig_brain;

segStruct.vol1 = orig;
clear orig;

segStruct.brainmask = segStruct.vol2;
clear orig_brain_mask;
segStruct.icv = length(find((segStruct.topLabel.imgs > 0) & (segStruct.topLabel.imgs~=22)))/1000;

segStruct.M_atl_to_vols = eye(4,4);

segStruct.anatomical_atlas.NumClasses = 23;
segStruct.feature_atlas.NumClasses = 23;
segStruct.numSpectra = 1;
segStruct.pvFactor = -1;
segStruct.pvList = zeros(0,5);

if DO_SAVE

    fprintf(fp, 'Saving %s/segStruct%s.mat ... ', pwd, strID);
    tic;
    save(['segStruct' strID '.mat'], 'segStruct');
    t = toc;
    totTime = totTime + t;
    fprintf(fp, 'done (%.1f s)\n', t);

else

    fprintf(fp, 'SKIPPING saving to segStruct.mat ... \n');

end
if DO_WRITEDICOM

    fprintf(fp, 'Saving DICOMs ... ', pwd, strID);
    tic;
    segdcmdir = [imgDir '//segdcm'];
    write_segdicom(segStruct.dcmhdr, segStruct,getColorMapCTX('gen'), segdcmdir, true, true, true, true, true);
    t = toc;
    totTime = totTime + t;
    fprintf(fp, 'done (%.1f s)\n', t);

else

    fprintf(fp, 'SKIPPING writing dicom files ... \n');

end

if DO_WRITESCDCM | DO_WRITEJPEG

    fprintf(fp, 'Saving JPEGS or SC DICOM ... ', pwd, strID);
    tic;
    write_seg2SCJPG(segStruct.dcmhdr, segStruct, getColorMapCTX('gen'), pwd, DO_WRITESCDCM, DO_WRITEJPEG);
    t = toc;
    totTime = totTime + t;
    fprintf(fp, 'done (%.1f s)\n', t);

else

    fprintf(fp, 'SKIPPING writing JPEGS/SC DCMS files ... \n');

end


if DO_SAVE
results = [segStruct.MI;segStruct.icv;segStruct.Left_cm3;segStruct.Right_cm3];
save results.txt results -ascii;
end

fprintf(fp, '---------------------------\n');
fprintf(fp, 'Total time: %.1f minutes\n', totTime/60);
fprintf(fp, '---------------------------\n');
fprintf(fp, 'JOBDONE\n');

if DISPLAY

    fprintf('\n');
    showSeg(segStruct);
    %keyboard

end

return

