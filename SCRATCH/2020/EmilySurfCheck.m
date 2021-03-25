clear; clc;

%% Amnestic MCI Amnestic VERSUS Controls
dir_loc = '/home/mmilmcdRSI/MetaData/ADNI/SurfGroupAvgs/AmnesticVNC';
fle_nme = 'Cortthickness_diff_tval_Amnestic_Normal Control-sm2819';

plt_dat{1} = fs_load_mgh( [ dir_loc '/' fle_nme '-' 'lh.mgh' ] );
plt_dat{2} = fs_load_mgh( [ dir_loc '/' fle_nme '-' 'rh.mgh' ] );

%
pcfg = [];

pcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/';
pcfg.dta_dir = [];
pcfg.out_pre_fix = 'EmTest_AmnesticMCI';

pcfg.plt_dta = { plt_dat{1} plt_dat{2} };

pcfg.fmr_col_map = {'neon blue' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
pcfg.fmr_rng_num = [];

pcfg.deg_fre = 350;

mmil_anat_surf_plot(pcfg)

%% Amnestic MCI Dysexecutive VERSUS Controls
dir_loc = '/home/mmilmcdRSI/MetaData/ADNI/SurfGroupAvgs/DysexecutiveVNC';
fle_nme = 'Cortthickness_diff_tval_Dysexecutive_Normal Control-sm2819';

plt_dat{1} = fs_load_mgh( [ dir_loc '/' fle_nme '-' 'lh.mgh' ] );
plt_dat{2} = fs_load_mgh( [ dir_loc '/' fle_nme '-' 'rh.mgh' ] );

%
pcfg = [];

pcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/';
pcfg.dta_dir = [];
pcfg.out_pre_fix = 'EmTest_DysexecutiveMCI';

pcfg.plt_dta = { plt_dat{1} plt_dat{2} };

pcfg.fmr_col_map = {'neon blue' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
pcfg.fmr_rng_num = [];

pcfg.deg_fre = 350;

mmil_anat_surf_plot(pcfg)

%% Pseudo-Don Amnestic MCI Dysexecutive VERSUS Controls
dir_loc = '/home/ekaestne/PROJECTS/OUTPUT/EmTest/SurfGroupAvgs/';
fle_nme = 'CortThick_pc_diff_tval_Amnestic_Normal Control-sm2819';

plt_dat{1} = fs_load_mgh( [ dir_loc '/' fle_nme '-' 'lh.mgh' ] );
plt_dat{2} = fs_load_mgh( [ dir_loc '/' fle_nme '-' 'rh.mgh' ] );

%
pcfg = [];

pcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/EmTest';
pcfg.dta_dir = [];
pcfg.out_pre_fix = 'EmTest_AmnesticMCI_PseudoDon_v2';

pcfg.plt_dta = { plt_dat{1} plt_dat{2} };

pcfg.fmr_col_map = {'neon blue' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
pcfg.fmr_rng_num = [];

pcfg.deg_fre = 460;

mmil_anat_surf_plot(pcfg)

%%
grp_lng = mmil_readtext('/home/ekaestne/PROJECTS/SUBJECTS/tmp/epd_age/ADNI_VisitInfo_cleaned.csv');
grp_lng = grp_lng(:,[ 1 2 5 ]);

analysis_outdir = 'analysis';
hemilist = {'lh','rh'};
sphsmooth = 2819;

dir_fsr = '/space/syn02/1/data/MMILDB/carrierm/fsurf';
outdir = sprintf( [ '/home/ekaestne/PROJECTS/OUTPUT/' '/' 'EmTest' '/' 'SurfGroupAvgs' ] );

dirlist_fs = dir( [ dir_fsr '/' 'FSURF*' ] );
[success,message]=mkdir(outdir);
fname = sprintf('%s/SurfGroupAvg_CortThick.log',outdir);
flog = fopen(fname,'wt');

Diags{1} = 'Amnestic';
Diags{2} = 'Normal Control';

Diagsbak = grp_lng(:,3);

total_sum = 0;
fprintf('%s: loading data...\n',mfilename);
fprintf(flog,'Studies included in surface group average:\n');
for d=1:length(Diags) % loop thru each diag    
  diag = Diags{d};
  S_ind = find(strcmp(Diagsbak,diag));
  for h = 1:length(hemilist) % loop thru each hemi
    hemi = hemilist{h};
    fprintf(flog,'diagnosis %s, hemi %s\n',diag,hemi);
    info{d}.info{h}.nsum = 0;
    info{d}.info{h}.vecsum = 0;
    info{d}.info{h}.vecsum2 = 0;
    info{d}.info{h}.data = nan(numel(S_ind),163842);
    for s = 1:length(S_ind) % loop thru each subject of specific diag
      SubjID = grp_lng{S_ind(s),2};
      n = regexp({dirlist_fs.name},['FSURF_' SubjID '_\d{8}.+_1$'],'match');
      n = [n{:}];
      if isempty(n)
          fprintf('%s: WARNING: no Recon container found for subject %s\n',mfilename,SubjID);
          continue;
      end
      ContainerDir = n{1};
      ContainerPath = sprintf('%s/%s',dir_fsr,ContainerDir);

      fname = sprintf('%s/%s/thickness-sphere-sm%d-%s.mgz',...
          ContainerPath,analysis_outdir,sphsmooth,hemi);
      if ~exist(fname,'file')
          fprintf('%s: WARNING: %s not found\n',mfilename,fname);
          continue;
      end;
      vec = mmil_rowvec(fs_load_mgh(fname));
      if isempty(vec)
          fprintf('%s: WARNING: unable to correctly read %s\n',mfilename,fname);
          continue;
      end
      info{d}.info{h}.data(s,:) = vec;
      info{d}.info{h}.nsum = info{d}.info{h}.nsum + 1;
      info{d}.info{h}.vecsum = info{d}.info{h}.vecsum + vec;
      info{d}.info{h}.vecsum2 = info{d}.info{h}.vecsum2 + vec.^2;
      fprintf(flog,'%s\n',ContainerDir);
    end
    nsum = info{d}.info{h}.nsum;
    fprintf(flog,'N for diagnosis %s, hemi %s: %d\n',diag,hemi,nsum);
    if nsum<=1
        fprintf('%s: WARNING: <=1 subjects included for diag %s, hemi %s\n',...
            mfilename,diag,hemi);
        nsum=2;
    end;
    info{d}.info{h}.vecmean = info{d}.info{h}.vecsum / (eps+nsum);
    info{d}.info{h}.vecstd = sqrt((nsum*info{d}.info{h}.vecsum2 - info{d}.info{h}.vecsum.^2)./(eps+nsum*(nsum-1)));
    info{d}.info{h}.vecstderr = info{d}.info{h}.vecstd/sqrt(eps+nsum);
  end
  total_sum = total_sum + nsum;
end
fprintf('%s: TOTAL n = %d\n',mfilename,total_sum);
fclose(flog);

save( '/home/ekaestne/PROJECTS/OUTPUT/EmTest/data_hold.mat' , 'info' )

fprintf('%s: saving group averages...\n',mfilename);
for d=1:length(Diags)
  for h = 1:length(hemilist)
    hemi = hemilist{h};
    vecmean = info{d}.info{h}.vecmean;
    vecstderr = info{d}.info{h}.vecstderr;
    vectval = vecmean./vecstderr;
    fname = sprintf('%s/CortThick_mean_%s-sm%d-%s.mgh',...
        outdir,Diags{d},sphsmooth,hemi);
    fs_save_mgh(vecmean,fname,eye(4));
    fname = sprintf('%s/CortThick_tval_%s-sm%d-%s.mgh',...
        outdir,Diags{d},sphsmooth,hemi);
    fs_save_mgh(vectval,fname,eye(4));
    fname = sprintf('%s/CortThick_stderr_%s-sm%d-%s.mgh',...
        outdir,Diags{d},sphsmooth,hemi);
    fs_save_mgh(vecstderr,fname,eye(4));
  end
end

contrastlist = [ 1 2 ];
for c = 1:size(contrastlist,1)
  i1 = contrastlist(c,1);
  i2 = contrastlist(c,2);
  for h = 1:length(hemilist)
    hemi = hemilist{h};

    % absolute difference
    vecdiff = info{i1}.info{h}.vecmean - info{i2}.info{h}.vecmean;
    vecstderr = sqrt(info{i1}.info{h}.vecstderr.^2+info{i2}.info{h}.vecstderr.^2);
    vectval = vecdiff./vecstderr;
    fname = sprintf('%s/CortThick_diff_mean_%s_%s-sm%d-%s.mgh',...
        outdir,Diags{i1},Diags{i2},sphsmooth,hemi);
    fs_save_mgh(vecdiff,fname,eye(4));
    fname = sprintf('%s/CortThick_diff_tval_%s_%s-sm%d-%s.mgh',...
        outdir,Diags{i1},Diags{i2},sphsmooth,hemi);
    fs_save_mgh(vectval,fname,eye(4));
    fname = sprintf('%s/CortThick_diff_stderr_%s_%s-sm%d-%s.mgh',...
        outdir,Diags{i1},Diags{i2},sphsmooth,hemi);
    fs_save_mgh(vecstderr,fname,eye(4));

    % fractional difference
    vecdiff = (info{i1}.info{h}.vecmean - info{i2}.info{h}.vecmean)./info{1}.info{h}.vecmean; % Difference in percent of NL
    vecstderr = sqrt((info{i1}.info{h}.vecstderr./info{1}.info{h}.vecmean).^2+...
        (info{i2}.info{h}.vecstderr./info{1}.info{h}.vecmean).^2);
    vectval = vecdiff./vecstderr;
    fname = sprintf('%s/CortThick_pc_diff_mean_%s_%s-sm%d-%s.mgh',...
        outdir,Diags{i1},Diags{i2},sphsmooth,hemi);
    fs_save_mgh(vecdiff,fname,eye(4));
    fname = sprintf('%s/CortThick_pc_diff_tval_%s_%s-sm%d-%s.mgh',...
        outdir,Diags{i1},Diags{i2},sphsmooth,hemi);
    fs_save_mgh(vectval,fname,eye(4));
    fname = sprintf('%s/CortThick_pc_diff_stderr_%s_%s-sm%d-%s.mgh',...
        outdir,Diags{i1},Diags{i2},sphsmooth,hemi);
    fs_save_mgh(vecstderr,fname,eye(4));
  end
end


