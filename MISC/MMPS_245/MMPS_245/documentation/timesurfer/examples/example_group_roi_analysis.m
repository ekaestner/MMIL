% group_roi_analysis

plotsubjflag = 1;
plotgroupflag = 1;
calcgroupflag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create input parameters

subjdir = '/home/dhagler/data/tp_roi_analysis/subjects';
labeldir = '/home/dhagler/data/tp_roi_analysis/labels/dhagler';
stcdir = '/home/dhagler/data/tp_roi_analysis/stcfiles';

% exclude subject nathanf because stc file has only 451 time points (sample=3.33)
% exclude subject clarka because waveforms are abnormal

subjnames = {...
'ben2'...
'christopher2'...
'josh2'...
'stephen2'...
'stewart2'...
'tim2'...
};
subjabrevs = {...
'bj'...
'ck'...
'jc'...
'sl'...
'sr'...
'ts'...
};

roinames = {...
'frontal'...
'SPL'...
};

stcprefixes = {...
'tpa_left'...
'tpa_right'...
'tpleft'...
'tpright'...
'tp-tpa_left'...
'tp-tpa_right'...
};
condnames = {...
'tpa left'...
'tpa right'...
'tp left'...
'tp right'...
'tp-tpa left'...
'tp-tpa right'...
};


contrasts = {...
[1 2]...
[3 4]...
[5 6]...
[3 1]...
[4 2]...
};

contrasts = {...
[3 4]...
[5 6]...
};

contrasts = {...
[3 4]...
};

surfname = 'orig';
hemilist = {'lh','rh'};

plottimes = [-100,400];
wformlimits = [0,6];
pvallimits = [1,5];

nsubs = length(subjnames);
nconds = length(stcprefixes);
ncontrasts = length(contrasts);
nhemis = length(hemilist);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

groupflag = 1;
tpoints = [];
nrois = [];
for h=1:nhemis
  hemi = hemilist{h};
  for c=1:nconds
    stcprefix = stcprefixes{c};
    condname = condnames{c};

    % create stcfile list
    stcfiles = [];
    for i=1:nsubs
      fullstcdir = sprintf('%s/%s_tp/sol',stcdir,subjabrevs{i});
      dirlist = dir(fullstcdir);
      stcfile = [];
      for j=1:length(dirlist)
        if ~dirlist(j).isdir &...
                       findstr(dirlist(j).name,stcprefix) &...
                       findstr(dirlist(j).name,hemi) &...
                       findstr(dirlist(j).name,'.stc')
           stcfile = sprintf('%s/%s',fullstcdir,dirlist(j).name);
        end;
      end;
      stcfiles{i} = stcfile;
      if isempty(stcfile)
        fprintf('%s: warning, stcfile not found for subj %s, prefix %s\n',...
          mfilename,subjnames{i},stcprefix);
        fprintf('%s:   looked for stc files in %s\n',...
          mfilename,fullstcdir);
      end;
    end;

    % create label list
    labelfiles = [];
    for i=1:nsubs
      for j=1:length(roinames)
        labelfiles(i).fnames{j} = sprintf('%s/%s/%s-%s.label',...
          labeldir,subjnames{i},hemi,roinames{j});
      end;
    end;

    % extract waveforms for ROIs
    matfile = sprintf('./subjwforms-%s-%s.mat',stcprefix,hemi);
    if ~exist(matfile,'file')
      fprintf('%s: extracting waveforms for cond %s...\n',...
        mfilename,stcprefix);
      results = ts_group_stc2roi(subjnames,stcfiles,labelfiles,...
        hemi,surfname,subjdir);
      if ~isempty(results)
        save(matfile,'results');
      end;
    else
      load(matfile);
    end;

    if isempty(tpoints)
      [tpoints,nrois] = size(results(1).wforms);
    else
      [tmp_tpoints,tmp_nrois] = size(results(1).wforms);
      if tmp_tpoints ~= tpoints
        fprintf('%s: warning: mismatch in number of timepoints for %s\n',...
          mfilename,stcprefix);
        groupflag = 0;
      end;
    end;

    % calculate average waveform
    matfile = sprintf('./avgwforms-%s-%s.mat',stcprefix,hemi);
    if ~exist(matfile,'file')
      fprintf('%s: calculating average waveform...\n',mfilename);
      time = results(1).time;
      avg = zeros(size(results(1).wforms));
      for i=1:nsubs
        [tmp_tpoints,tmp_nrois]=size(results(i).wforms);
        if tmp_tpoints ~= tpoints
          fprintf('%s: error: mismatch in number of timepoints for subj %s\n',...
            mfilename,subjnames{i});
          return;
        end;
        avg = avg + results(i).wforms;
      end;
      avg = avg/length(subjnames);
      save(matfile,'avg','time');
    else
      load(matfile);
    end;

    if plotsubjflag & 0
      fprintf('%s: plotting subject waveforms...\n',mfilename);
      clf;
      for i=1:nsubs
        subplot(nsubs,1,i);
        plot(results(i).time,results(i).wforms);
        set(gca,'YLim',wformlimits);
        set(gca,'XLim',plottimes);
        title(sprintf('subject %s',subjnames{i}));
      end;
      fname = sprintf('subjwforms-%s-%s.jpg',stcprefix,hemi);
      print('-djpeg',fname);

%      fname = sprintf('subjwforms-%s-%s.eps',stcprefix,hemi);
%      print('-deps',fname);

%      fname = sprintf('subjwforms-%s-%s.ai',stcprefix,hemi);
%      print('-dill',fname);
    end;
    
    if plotgroupflag & 0
      fprintf('%s: plotting group avg waveforms...\n',mfilename);
      clf;
      plot(time,avg);
      set(gca,'YLim',wformlimits);
      set(gca,'XLim',plottimes);
      title(sprintf('%s group avg',condname));
      legend(roinames,'Location','East');
      fname = sprintf('avgwforms-%s-%s.jpg',stcprefix,hemi);
      print('-djpeg',fname);
    end;
  end;
end;

% compile results
if groupflag
  fprintf('%s: compiling results...\n',mfilename);
  all_results = [];
  for h=1:nhemis
    hemi = hemilist{h};
    all_results{h}.data = zeros(nconds,tpoints,nrois,nsubs);
    for c=1:nconds
      stcprefix = stcprefixes{c};
      % load subject waveforms
      matfile = sprintf('./subjwforms-%s-%s.mat',stcprefix,hemi);
      results = [];
      load(matfile);
      time = results(1).time;
      for s=1:nsubs
        all_results{h}.data(c,:,:,s) = results(s).wforms;
      end;
    end;
  end;
else
  fprintf('%s: unable to do group stats (mismatch in tpoints)\n',...
    mfilename);
  return;
end;

if plotsubjflag
  fprintf('%s: plotting group avg waveforms...\n',mfilename);
  for s=1:nsubs
    subjname = subjnames{s};
    for r=1:nrois
      roiname = roinames{r};
      for i=1:ncontrasts
        tmp_cont = contrasts{i};
        cond1 = tmp_cont(1);
        cond2 = tmp_cont(2);
        clf;
        for h=1:nhemis
          hemi = hemilist{h};
          tmp_data = squeeze(all_results{h}.data(tmp_cont,:,r,s));
          subplot(1,2,h);
          plot(time,tmp_data);
          set(gca,'YLim',wformlimits);
          set(gca,'XLim',plottimes);
          title(sprintf('%s %s %s',hemi,roiname,subjname));
          legend(condnames{tmp_cont},'Location','NorthEast');
        end;
        fname = sprintf('subjwforms-%s-%sVS%s-%s.jpg',subjname,stcprefixes{tmp_cont},roiname);
        print('-djpeg',fname);
      end;
    end;
  end;
end;


if plotgroupflag
  fprintf('%s: plotting group avg waveforms...\n',mfilename);
  for r=1:nrois
    roiname = roinames{r};
    for i=1:ncontrasts
      tmp_cont = contrasts{i};
      cond1 = tmp_cont(1);
      cond2 = tmp_cont(2);
      clf;
      for h=1:nhemis
        hemi = hemilist{h};
        tmp_data = squeeze(all_results{h}.data(tmp_cont,:,r,:));
        tmp_mean = mean(tmp_data,3);
        subplot(1,2,h);
        plot(time,tmp_mean);
        set(gca,'YLim',wformlimits);
        set(gca,'XLim',plottimes);
        title(sprintf('%s %s group avg',hemi,roiname));
        legend(condnames{tmp_cont},'Location','NorthEast');
      end;
      fname = sprintf('avgwforms-%sVS%s-%s.jpg',stcprefixes{tmp_cont},roiname);
      print('-djpeg',fname);
    end;
  end;
end;

if ~calcgroupflag, return; end;

%% one-sample ttests
matfile = sprintf('./ttests.mat');
if ~exist(matfile,'file')
  fprintf('%s: calculating one-sample t-tests...\n',mfilename);
  for h=1:nhemis
    hemi = hemilist{h};
    all_results{h}.mean = mean(all_results{h}.data,4);
    all_results{h}.stdv = std(all_results{h}.data,0,4);
    all_results{h}.stat = all_results{h}.mean./(all_results{h}.stdv/sqrt(nsubs));
    all_results{h}.pval = 1-tcdf(all_results{h}.stat,nsubs-1);
  end;
  save(matfile,'all_results');
else
  load(matfile);
end;

%% two-sample ttests
matfile = sprintf('./twosamp_ttests.mat');
if ~exist(matfile,'file')
  fprintf('%s: calculating 2-sample paired ttests...\n',mfilename);
  for h=1:nhemis
    twosamptest = [];
    for i=1:ncontrasts
      tmp_cont = contrasts{i};
      cond1 = tmp_cont(1);
      cond2 = tmp_cont(2);
      data1 = squeeze(all_results{h}.data(cond1,:,:,:));
      data2 = squeeze(all_results{h}.data(cond2,:,:,:));
      diff = data1 - data2;
      diff_mean = mean(diff,3);
      diff_std = std(diff,0,3);
      diff_stat = diff_mean./(diff_std/sqrt(nsubs));
      diff_pval = 1-tcdf(diff_stat,nsubs-1);
      twosampttest(i).condname1 = condnames{cond1};
      twosampttest(i).condname2 = condnames{cond2};
      twosampttest(i).diff = diff;
      twosampttest(i).mean = diff_mean;
      twosampttest(i).std = diff_std;
      twosampttest(i).stat = diff_stat;
      twosampttest(i).pval = diff_pval;
    end;
    all_results{h}.twosampttest = twosampttest;
  end;
  save(matfile,'all_results');
else
  load(matfile);
end;

%% sign-rank tests
matfile = sprintf('./signrank.mat');
if ~exist(matfile,'file')
  fprintf('%s: calculating Wilcoxon signed rank tests...\n',mfilename);
  for h=1:nhemis
    signranktest = [];
    for i=1:ncontrasts
      tmp_cont = contrasts{i};
      cond1 = tmp_cont(1);
      cond2 = tmp_cont(2);
      data1 = squeeze(all_results{h}.data(cond1,:,:,:));
      data2 = squeeze(all_results{h}.data(cond2,:,:,:));
      diff = data1 - data2;
      diff_mean = mean(diff,3);
      diff_std = std(diff,0,3);
      diff_median = median(diff,3);
      diff_stat = zeros(tpoints,nrois);
      diff_pval = zeros(tpoints,nrois);
      for t=1:tpoints
        for r=1:nrois
          tmp_data = squeeze(diff(t,r,:));
          [pval,hyp,stats] = signrank(tmp_data);
          diff_stat(t,r) = stats.signedrank;
          diff_pval(t,r) = pval;
        end;
      end;
      signranktest(i).condname1 = condnames{cond1};
      signranktest(i).condname2 = condnames{cond2};
      signranktest(i).diff = diff;
      signranktest(i).median = diff_median;
      signranktest(i).mean = diff_mean;
      signranktest(i).std = diff_std;
      signranktest(i).stat = diff_stat;
      signranktest(i).pval = diff_pval;
    end;
    all_results{h}.signranktest = signranktest;
  end;
  save(matfile,'all_results');
else
  load(matfile);
end;

%% 3-way anova (hemi fixed effect, condition fixed effect, subject random effect)
matfile = sprintf('./anova.mat');
if ~exist(matfile,'file')
  fprintf('%s: calculating 3-way ANOVAs...\n',mfilename);
  anova3 = [];
  for i=1:ncontrasts
    tmp_cont = contrasts{i};
    anova3(i).hemi_pval = zeros(tpoints,nrois);
    anova3(i).cond_pval = zeros(tpoints,nrois);
    anova3(i).subj_pval = zeros(tpoints,nrois);
    anova3(i).hemiXcond_pval = zeros(tpoints,nrois);
    anova3(i).hemiXsubject_pval = zeros(tpoints,nrois);
    anova3(i).condXsubject_pval = zeros(tpoints,nrois);
    anova3(i).condname1 = condnames{tmp_cont(1)};
    anova3(i).condname2 = condnames{tmp_cont(2)};
    for r=1:nrois
      roiname = roinames{r};
      fprintf('%s: roi %s\n',mfilename,roiname);
      for t=1:tpoints
        % construct vector of values
        nvals = nhemis*2*nsubs;
        Y = zeros(nvals,1);
        g_hemi = [];
        g_cond = [];
        g_subj = [];
        n=1;
        for h=1:nhemis
          for c=1:2
            cond = tmp_cont(c);
            for s=1:nsubs
              Y(n) = all_results{h}.data(cond,t,r,s);
              g_hemi(n) = h;
              g_cond(n) = cond;
              g_subj(n) = s;
              n=n+1;
            end;
          end;
        end;
        [p,table,stats]=anovan(Y,{g_hemi g_cond g_subj},...
          'model','full','random',[3],...
          'varnames',{'Hemi','Cond','Subject'},...
          'display','off');
        anova3(i).hemi_pval(t,r) = p(1);
        anova3(i).cond_pval(t,r) = p(2);
        anova3(i).subj_pval(t,r) = p(3);
        anova3(i).hemiXcond_pval(t,r) = p(4);
        anova3(i).hemiXsubject_pval(t,r) = p(5);
        anova3(i).condXsubject_pval(t,r) = p(6);
      end;
    end;
  end;
  save(matfile,'anova3');
else
  load(matfile);
end;


if plotgroupflag
  fprintf('%s: plotting t-test results...\n',mfilename);
  for h=1:nhemis
    hemi = hemilist{h};
    for c=1:nconds
      condname = condnames{c};
      clf;
      tmp_stat = squeeze(all_results{h}.stat(c,:,:));
      tmp_pval = -log10(squeeze(all_results{h}.pval(c,:,:)));
      tmp_mean = squeeze(all_results{h}.mean(c,:,:));
      subplot(3,1,1)
      plot(time,tmp_mean);
      title('mean');
      set(gca,'YLim',wformlimits);
      set(gca,'XLim',plottimes);
      subplot(3,1,2)
      plot(time,tmp_stat);
      title('stat');
      set(gca,'XLim',plottimes);
      subplot(3,1,3)
      plot(time,tmp_pval);
      title('-log10 pval');
      set(gca,'XLim',plottimes);
      set(gca,'Ylim',pvallimits);
      fname = sprintf('ttests-%s-%s.jpg',stcprefix,hemi);
      print('-djpeg',fname);
    end;
  end;
end;

if plotgroupflag
  fprintf('%s: plotting two sample paired t-test results...\n',mfilename);
  for h=1:nhemis
    hemi = hemilist{h};
    for i=1:ncontrasts
      clf;
      condname1 = all_results{h}.twosampttest(i).condname1;
      condname2 = all_results{h}.twosampttest(i).condname2;
      tmp_stat = all_results{h}.twosampttest(i).stat;
      tmp_pval = -log10(all_results{h}.twosampttest(i).pval);
      tmp_mean = all_results{h}.twosampttest(i).mean;
      subplot(3,1,1)
      plot(time,tmp_mean);
      title('mean');
      set(gca,'YLim',wformlimits);
      set(gca,'XLim',plottimes);
      subplot(3,1,2)
      plot(time,tmp_stat);
      title('stat');
      set(gca,'XLim',plottimes);
      subplot(3,1,3)
      plot(time,tmp_pval);
      title('-log10 pval');
      set(gca,'XLim',plottimes);
      set(gca,'Ylim',pvallimits);
      fname = sprintf('twosampttests-%sVS%s-%s.jpg',condname1,condname2,hemi);
      print('-djpeg',fname);
    end;
  end;
end;

if plotgroupflag
  fprintf('%s: plotting Wilcoxon signed rank test results...\n',mfilename);
  for h=1:nhemis
    hemi = hemilist{h};
    for i=1:ncontrasts
      clf;
      condname1 = all_results{h}.signranktest(i).condname1;
      condname2 = all_results{h}.signranktest(i).condname2;
      tmp_stat = all_results{h}.signranktest(i).stat;
      tmp_pval = -log10(all_results{h}.signranktest(i).pval);
      tmp_median = all_results{h}.signranktest(i).median;
      tmp_mean = all_results{h}.signranktest(i).mean;
      subplot(3,1,1)
      plot(time,tmp_median);
      title('median');
%      plot(time,tmp_mean);
%      title('mean');
      set(gca,'YLim',wformlimits);
      set(gca,'XLim',plottimes);
      subplot(3,1,2)
      plot(time,tmp_stat);
      title('stat');
      set(gca,'XLim',plottimes);
      subplot(3,1,3)
      plot(time,tmp_pval);
      title('-log10 pval');
      set(gca,'XLim',plottimes);
      set(gca,'Ylim',pvallimits);
      fname = sprintf('signrank-%sVS%s-%s.jpg',condname1,condname2,hemi);
      print('-djpeg',fname);
    end;
  end;
end;

if plotgroupflag
  fprintf('%s: plotting 3-way ANOVA results...\n',mfilename);
  for i=1:ncontrasts
    clf;
    condname1 = anova3(i).condname1;
    condname2 = anova3(i).condname2;

    tmp_hemi_pval = abs(anova3(i).hemi_pval);
    tmp_cond_pval = abs(anova3(i).cond_pval);
    tmp_subj_pval = abs(anova3(i).subj_pval);
    tmp_hemiXcond_pval = abs(anova3(i).hemiXcond_pval);
    tmp_hemiXsubject_pval = abs(anova3(i).hemiXsubject_pval);
    tmp_condXsubject_pval = abs(anova3(i).condXsubject_pval);

    tmp_hemi_pval(find(tmp_hemi_pval==0)) = 10^-10;
    tmp_cond_pval(find(tmp_cond_pval==0)) = 10^-10;
    tmp_subj_pval(find(tmp_subj_pval==0)) = 10^-10;
    tmp_hemiXcond_pval(find(tmp_hemiXcond_pval==0)) = 10^-10;
    tmp_hemiXsubject_pval(find(tmp_hemiXsubject_pval==0)) = 10^-10;
    tmp_condXsubject_pval(find(tmp_condXsubject_pval==0)) = 10^-10;

    tmp_hemi_pval = -log10(tmp_hemi_pval);
    tmp_cond_pval = -log10(tmp_cond_pval);
    tmp_subj_pval = -log10(tmp_subj_pval);
    tmp_hemiXcond_pval = -log10(tmp_hemiXcond_pval);
    tmp_hemiXsubject_pval = -log10(tmp_hemiXsubject_pval);
    tmp_condXsubject_pval = -log10(tmp_condXsubject_pval);

    subplot(6,1,1)
    plot(time,tmp_hemi_pval);
    title('-log10 hemi pval');
    set(gca,'XLim',plottimes);
    set(gca,'Ylim',pvallimits);

    subplot(6,1,2)
    plot(time,tmp_cond_pval);
    title('-log10 cond pval');
    set(gca,'XLim',plottimes);
    set(gca,'Ylim',pvallimits);

    subplot(6,1,3)
    plot(time,tmp_subj_pval);
    title('-log10 subj pval');
    set(gca,'XLim',plottimes);
    set(gca,'Ylim',pvallimits);

    subplot(6,1,4)
    plot(time,tmp_hemiXcond_pval);
    title('-log10 hemiXcond pval');
    set(gca,'XLim',plottimes);
    set(gca,'Ylim',pvallimits);

    subplot(6,1,5)
    plot(time,tmp_hemiXsubject_pval);
    title('-log10 hemiXsubject pval');
    set(gca,'XLim',plottimes);
    set(gca,'Ylim',pvallimits);

    subplot(6,1,6)
    plot(time,tmp_condXsubject_pval);
    title('-log10 condXsubject pval');
    set(gca,'XLim',plottimes);
    set(gca,'Ylim',pvallimits);

    fname = sprintf('avova3-%sVS%s.jpg',condname1,condname2);
    print('-djpeg',fname);
  end;
end;

%% todo: integrate time ranges, then do stats


