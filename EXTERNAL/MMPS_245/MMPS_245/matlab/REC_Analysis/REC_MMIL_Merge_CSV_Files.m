function REC_MMIL_Merge_CSV_Files(ProjID,normflag)
%function REC_MMIL_Merge_CSV_Files(ProjID,normflag)
%
% Required Inputs:
% ProjID - name of Recharge project to run 
% normflag - for PET, normalize to pons (default = 1)
%    (0 = absolute values, 1 = norm to pons, 2 = both)
%
%
% Early Mod: 06/12/09 by Alain Koyama
% Last Mod:  08/16/10 by Cooper Roddey
%

% todo: merge DTI data

if ~exist('ProjID','var'), ProjID = []; end;
if ~exist('normflag','var'), normflag = 1; end;

if isempty(ProjID)
   error('Empty ProjID');
end

[RootDirs,ProjInfo] = REC_MMIL_RootDirs(ProjID);
homedir = RootDirs.home;

for p=1:length(ProjInfo) % start project loop
  % if no PET, normflag = 0;
  if isempty(ProjInfo(p).TMP_PET_RootDir)
    normflag = 0; 
  end;
  switch normflag
    case 0, normflags = 0;
    case 1, normflags = 1;
    case 2, normflags = [0 1];
  end
  for nn=1:length(normflags)
    normflag = normflags(nn);
    ProjID = ProjInfo(p).ProjID;
    metadatarootdir = sprintf('%s/MetaData/%s',homedir,ProjInfo(p).ProjID);
    if normflag
      outmetadatafname = sprintf('%s/REC_%s_Complete_norm.csv',metadatarootdir,ProjID);
    else
      outmetadatafname = sprintf('%s/REC_%s_Complete.csv',metadatarootdir,ProjID);
    end
    fname = sprintf('%s/MetaData/%s/REC_%s_StudyInfo.csv',homedir,ProjID,ProjID);
    if exist(fname,'file') % load subject meta data if it exists
      ptdata = mmil_readtext(fname,',','','"');
    else
      fprintf('%s: ERROR: %s not found. Please generate StudyInfo using REC_MMIL_Import_StudyInfo\n',mfilename,fname);
      return;
    end
    SubjIDs = ptdata(:,1);
    SubjInfo = cell(size(ptdata,1),1); % all data goes here
    for i=1:size(ptdata,1) % copy info from studyinfo
      for j=1:3 % copy 1st 3 columns (subject ID, studydate, group)
        SubjInfo{i,j} = ptdata{i,j};
      end
    end
    h=size(SubjInfo,2);

    % read subcortical volumes
    metadatafname = sprintf('%s/MetaData/%s/ROI_Summaries/REC_%s_SubCortVolume.csv',...
      homedir,ProjID,ProjID);
    metadata = mmil_readtext(metadatafname,',','','"');
    h_start = h;
    n = regexpi(metadata(1,:),'StudyDate'); 
    n = find(~cellfun('isempty', n)==1); n = n(1)+1;
    for i=n:size(metadata,2) % copy SubCort headers
      h = h + 1;
      SubjInfo{1,h} = sprintf('SUBCORT_%s',metadata{1,i});
    end
    for i=2:size(metadata,1) % loop thru each subject in SubCort file
      h = h_start;
      SubjID = metadata{i,1};
      S_ind = find(strcmp(SubjIDs,SubjID));
      if isempty(S_ind)
        fprintf('%s: ERROR Subject %s not found in StudyInfo file %s\n',...
          mfilename, SubjID, fname);
        return;
      end
      for j=n:size(metadata,2) % loop thru each data column in SubCort file
        h = h + 1;
        SubjInfo{S_ind,h} = metadata{i,j};
      end
    end
    % read cortical volumes
    metadatafname = sprintf('%s/MetaData/%s/ROI_Summaries/REC_%s_CortThick.csv',...
      homedir,ProjID,ProjID);
    metadata = mmil_readtext(metadatafname,',','','"');
    h_start = h;
    n = regexpi(metadata(1,:),'StudyDate'); 
    n = find(~cellfun('isempty', n)==1); n = n(1)+1;
    for i=n:size(metadata,2) % copy SubCort headers
      h = h + 1;
      SubjInfo{1,h} = sprintf('CORT_%s',metadata{1,i});
    end
    for i=2:size(metadata,1) % loop thru each subject in SubCort file
      h = h_start;
      SubjID = metadata{i,1};
      S_ind = find(strcmp(SubjIDs,SubjID));
      for j=n:size(metadata,2) % loop thru each data column in SubCort file
        h = h + 1;
        SubjInfo{S_ind,h} = metadata{i,j};
      end
    end

    % read PET aseg
    if normflag
      metadatafnames = dir(sprintf('%s/MetaData/%s/ROI_Summaries/REC_%s_PET_norm_*aseg.csv',...
        homedir,ProjID,ProjID));
    else
      metadatafnames = dir(sprintf('%s/MetaData/%s/ROI_Summaries/REC_%s_PET*aseg.csv',...
        homedir,ProjID,ProjID));
      n = regexpi({metadatafnames.name},'norm'); % find series but exclude 'norm'
      n = ~cellfun('isempty', n);  n = ~n; n = find(n==1);
      metadatafnames = metadatafnames(n);
    end
    volnames = {};
    if isempty(metadatafnames) % find PET series for this project
%       fprintf('%s: Warning: No PET aseg files found\n',mfilename);
    else
      volnames = regexp({metadatafnames.name},['(?<=REC_' ProjID '_PET_).+(?=_aseg\.csv)'],'match');
    end

    for volnums=1:length(metadatafnames) % loop thru each PET series
      metadatafname = sprintf('%s/MetaData/%s/ROI_Summaries/%s',...
        homedir,ProjID,metadatafnames(volnums).name);
      metadata = mmil_readtext(metadatafname,',','','"');
      h_start = h;
      n = regexpi(metadata(1,:),'StudyDate'); 
      n = find(~cellfun('isempty', n)==1); n = n(1)+1;
      for i=n:size(metadata,2) % copy PET aseg headers
        h = h + 1;
        SubjInfo{1,h} = sprintf('PET_%s_%s',char(volnames{volnums}),metadata{1,i});
      end
      for i=2:size(metadata,1) % loop thru each subject in PET aseg file
        h = h_start;
        SubjID = metadata{i,1};
        S_ind = find(strcmp(SubjIDs,SubjID));
        for j=n:size(metadata,2) % loop thru each data column in PET aseg file
          h = h + 1;
          SubjInfo{S_ind,h} = metadata{i,j};
        end
      end
    end

    % read PET aparc
    if normflag
      metadatafnames = dir(sprintf('%s/MetaData/%s/ROI_Summaries/REC_%s_PET_norm_*aparc.csv',...
        homedir,ProjID,ProjID));
    else
      metadatafnames = dir(sprintf('%s/MetaData/%s/ROI_Summaries/REC_%s_PET*aparc.csv',...
        homedir,ProjID,ProjID));
      n = regexpi({metadatafnames.name},'norm'); % find series but exclude 'norm'
      n = ~cellfun('isempty', n);  n = ~n; n = find(n==1);
      metadatafnames = metadatafnames(n);
    end
    volnames = {};
    if isempty(metadatafnames) % find PET series for this project
%       fprintf('%s: Warning: No PET aparc files found\n',mfilename);
    else
      volnames = regexp({metadatafnames.name},['(?<=REC_' ProjID '_PET_).+(?=_aparc\.csv)'],'match');
    end

    for volnums=1:length(metadatafnames) % loop thru each PET series
      metadatafname = sprintf('%s/MetaData/%s/ROI_Summaries/%s',...
        homedir,ProjID,metadatafnames(volnums).name);
      metadata = mmil_readtext(metadatafname,',','','"');
      h_start = h;
      n = regexpi(metadata(1,:),'StudyDate'); 
      n = find(~cellfun('isempty', n)==1); n = n(1)+1;
      for i=n:size(metadata,2) % copy PET aseg  headers
        h = h + 1;
        SubjInfo{1,h} = sprintf('PET_%s_%s',char(volnames{volnums}),metadata{1,i});
      end
      for i=2:size(metadata,1) % loop thru each subject in PET aseg file
        h = h_start;
        SubjID = metadata{i,1};
        S_ind = find(strcmp(SubjIDs,SubjID));
        for j=n:size(metadata,2) % loop thru each data column in PET aseg file
          h = h + 1;
          SubjInfo{S_ind,h} = metadata{i,j};
        end
      end
    end
    
    % add all patient metadata
    ind = size(SubjInfo,2); % col # to start pasting data
    for i=2:size(ptdata,2) % copy metadata headers
      SubjInfo{1,i+ind-1} = ptdata{1,i};
    end

    for i=2:size(ptdata,1) % loop thru each subject in metadata file
      for j=2:size(ptdata,2) % loop thru each data column in metadata file
        SubjInfo{i,j+ind-1} = ptdata{i,j};
      end
    end

    for i=1:size(SubjInfo,1) % replace all empty cells to .
      for j=1:size(SubjInfo,2)
        if isempty(SubjInfo{i,j}), SubjInfo{i,j} = '.'; end;
      end
    end

    write_cellarray_csv(outmetadatafname,SubjInfo);

    %%%%%% calculate averages %%%%%%

    SubjInfoAvg = {};
    for i=1:size(ptdata,1) % copy info from studyinfo
      for j=1:3 % copy 1st 3 columns
        SubjInfoAvg{i,j} = ptdata{i,j};
      end
    end
    if normflag
      outmetadatafnameavg = sprintf('%s/REC_%s_Complete_norm_Avg.csv',metadatarootdir,ProjID);
    else
      outmetadatafnameavg = sprintf('%s/REC_%s_Complete_Avg.csv',metadatarootdir,ProjID);
    end


    % subcort data to be averaged
    metadatafname = sprintf('%s/MetaData/%s/ROI_Summaries/REC_%s_SubCortVolume.csv',...
      homedir,ProjID,ProjID);
    metadata = mmil_readtext(metadatafname,',','','"');

    h=size(SubjInfoAvg,2);
    n = regexpi(metadata(1,:),'StudyDate'); 
    n = find(~cellfun('isempty', n)==1); n = n(1)+1;
    for j=n:size(metadata,2) % loop thru data
      idx_left = regexp(metadata{1,j},'Left','once');
      idx_right = regexp(metadata{1,j},'Right','once');
      if ~isempty(idx_left) % convert LH data col to Avg
        h=h+1;
        newheader = sprintf('SUBCORT_Avg%s',metadata{1,j}(idx_left+4:end));
        SubjInfoAvg{1,h} = newheader;
        % loop through all subjects and avg values for this col
        for k=2:size(metadata,1)
          rightheader = sprintf('%sRight%s',metadata{1,j}(1:idx_left-1),metadata{1,j}(idx_left+4:end));
          headers = {metadata{1,j:end}}; % remaining subcort headers
          idx=find(strcmp(headers,rightheader)); % find corresponding col for RH data
          lh_val = metadata{k,j};
          rh_val = metadata{k,j+idx-1};
          if ~isnumeric(lh_val) | ~isnumeric(rh_val)
            SubjInfoAvg{k,h} = '';
          else
            SubjInfoAvg{k,h} = (lh_val + rh_val)/2; % average LH and RH value
          end
        end
        continue;
      elseif ~isempty(idx_right) % omit RH data cols
        continue;
      else
        h=h+1;
        newheader = sprintf('SUBCORT_%s',metadata{1,j});
        SubjInfoAvg{1,h} = newheader;
        % loop through all subjects and copy non-avg values for this col
        for k=2:size(metadata,1)
          SubjInfoAvg{k,h} = metadata{k,j};
        end
      end
    end

    % cort data to be averaged
    metadatafname = sprintf('%s/MetaData/%s/ROI_Summaries/REC_%s_CortThick.csv',...
      homedir,ProjID,ProjID);
    metadata = mmil_readtext(metadatafname,',','','"');
    n = regexpi(metadata(1,:),'StudyDate'); 
    n = find(~cellfun('isempty', n)==1); n = n(1)+1;
    for j=n:size(metadata,2) % loop thru data
      idx_left = regexp(metadata{1,j},'lh','once');
      idx_right = regexp(metadata{1,j},'rh','once');
      if ~isempty(idx_left) % convert LH data col to Avg
        h=h+1;
        newheader = sprintf('CORT_Avg%s',metadata{1,j}(idx_left+2:end));
        SubjInfoAvg{1,h} = newheader;
        % loop through all subjects and avg values for this col
        for k=2:size(metadata,1)
          rightheader = sprintf('%srh%s',metadata{1,j}(1:idx_left-1),metadata{1,j}(idx_left+2:end));
          headers = {metadata{1,j:end}}; % remaining cort headers
          idx=find(strcmp(headers,rightheader)); % find corresponding col for RH data
          lh_val = metadata{k,j};
          rh_val = metadata{k,j+idx-1};
          if ~isnumeric(lh_val) | ~isnumeric(rh_val)
            SubjInfoAvg{k,h} = '';
          else
            SubjInfoAvg{k,h} = (lh_val + rh_val)/2; % average LH and RH value
          end
        end
        continue;
      elseif ~isempty(idx_right) % omit RH data cols
        continue;
      else
        h=h+1;
        newheader = sprintf('CORT_%s',metadata{1,j});
        SubjInfoAvg{1,h} = newheader;
        % loop through all subjects and copy non-avg values for this col
        for k=2:size(metadata,1)
          SubjInfoAvg{k,h} = metadata{k,j};
        end
      end
    end

    % PET aseg data to be averaged
    if normflag
      metadatafnames = dir(sprintf('%s/MetaData/%s/ROI_Summaries/REC_%s_PET_norm_*aseg.csv',...
        homedir,ProjID,ProjID));
    else
      metadatafnames = dir(sprintf('%s/MetaData/%s/ROI_Summaries/REC_%s_PET*aseg.csv',...
        homedir,ProjID,ProjID));
      n = regexpi({metadatafnames.name},'norm'); % find series but exclude 'norm'
      n = ~cellfun('isempty', n);  n = ~n; n = find(n==1);
      metadatafnames = metadatafnames(n);
    end
    volnames = {};
    if isempty(metadatafnames) % find PET series for this project
%       fprintf('%s: Warning: No PET aseg files found\n',mfilename);
    else
      volnames = regexp({metadatafnames.name},['(?<=REC_' ProjID '_PET_).+(?=_aseg\.csv)'],'match');
    end

    for volnums=1:length(metadatafnames) % loop thru each PET series
      metadatafname = sprintf('%s/MetaData/%s/ROI_Summaries/%s',...
        homedir,ProjID,metadatafnames(volnums).name);
      metadata = mmil_readtext(metadatafname,',','','"');
      n = regexpi(metadata(1,:),'StudyDate'); 
      n = find(~cellfun('isempty', n)==1); n = n(1)+1;
      for j=n:size(metadata,2) % loop thru data
        idx_left = regexp(metadata{1,j},'Left','once');
        idx_right = regexp(metadata{1,j},'Right','once');
        if ~isempty(idx_left) % convert LH data col to Avg
          h=h+1;
          newheader = sprintf('PET_%s_Avg%s',char(volnames{volnums}),metadata{1,j}(idx_left+4:end));
          SubjInfoAvg{1,h} = newheader;
          % loop through all subjects and avg values for this col
          for k=2:size(metadata,1)
            rightheader = sprintf('%sRight%s',metadata{1,j}(1:idx_left-1),metadata{1,j}(idx_left+4:end));
            headers = {metadata{1,j:end}}; % remaining headers
            idx=find(strcmp(headers,rightheader)); % find corresponding col for RH data
            lh_val = metadata{k,j};
            rh_val = metadata{k,j+idx-1};
            if ~isnumeric(lh_val) | ~isnumeric(rh_val)
              SubjInfoAvg{k,h} = '';
            else
              SubjInfoAvg{k,h} = (lh_val + rh_val)/2; % average LH and RH value
            end
          end
          continue;
        elseif ~isempty(idx_right) % omit RH data cols
          continue;
        else
          h=h+1;
          newheader = sprintf('PET_%s_%s',char(volnames{volnums}),metadata{1,j});
          SubjInfoAvg{1,h} = newheader;
          % loop through all subjects and copy non-avg values for this col
          for k=2:size(metadata,1)
            SubjInfoAvg{k,h} = metadata{k,j};
          end
        end
      end
    end

    % PET aparc data to be averaged
    if normflag
      metadatafnames = dir(sprintf('%s/MetaData/%s/ROI_Summaries/REC_%s_PET_norm_*aparc.csv',...
        homedir,ProjID,ProjID));
    else
      metadatafnames = dir(sprintf('%s/MetaData/%s/ROI_Summaries/REC_%s_PET*aparc.csv',...
        homedir,ProjID,ProjID));
      n = regexpi({metadatafnames.name},'norm'); % find series but exclude 'norm'
      n = ~cellfun('isempty', n);  n = ~n; n = find(n==1);
      metadatafnames = metadatafnames(n);
    end
    volnames = {};
    if isempty(metadatafnames) % find PET series for this project
%       fprintf('%s: Warning: No PET aparc files found\n',mfilename);
    else
      volnames = regexp({metadatafnames.name},['(?<=REC_' ProjID '_PET_).+(?=_aparc\.csv)'],'match');
    end

    for volnums=1:length(metadatafnames) % loop thru each PET series
      metadatafname = sprintf('%s/MetaData/%s/ROI_Summaries/%s',...
        homedir,ProjID,metadatafnames(volnums).name);
      metadata = mmil_readtext(metadatafname,',','','"');
      n = regexpi(metadata(1,:),'StudyDate'); 
      n = find(~cellfun('isempty', n)==1); n = n(1)+1;
      for j=n:size(metadata,2) % loop thru data
        idx_left = regexp(metadata{1,j},'lh','once');
        idx_right = regexp(metadata{1,j},'rh','once');
        if ~isempty(idx_left) % convert LH data col to Avg
          h=h+1;
          newheader = sprintf('PET_%s_Avg%s',char(volnames{volnums}),metadata{1,j}(idx_left+2:end));
          SubjInfoAvg{1,h} = newheader;
          % loop through all subjects and avg values for this col
          for k=2:size(metadata,1)
            rightheader = sprintf('%srh%s',metadata{1,j}(1:idx_left-1),metadata{1,j}(idx_left+2:end));
            headers = {metadata{1,j:end}}; % remaining headers
            idx=find(strcmp(headers,rightheader)); % find corresponding col for RH data
            lh_val = metadata{k,j};
            rh_val = metadata{k,j+idx-1};
            if ~isnumeric(lh_val) | ~isnumeric(rh_val)
              SubjInfoAvg{k,h} = '';
            else
              SubjInfoAvg{k,h} = (lh_val + rh_val)/2; % average LH and RH value
            end
          end
          continue;
        elseif ~isempty(idx_right) % omit RH data cols
          continue;
        else
          h=h+1;
          newheader = sprintf('PET_%s_%s',char(volnames{volnums}),metadata{1,j});
          SubjInfoAvg{1,h} = newheader;
          % loop through all subjects and copy non-avg values for this col
          for k=2:size(metadata,1)
            SubjInfoAvg{k,h} = metadata{k,j};
          end
        end
      end
    end

    % add all patient metadata
    ind = size(SubjInfoAvg,2); % col # to start pasting data
    for i=2:size(ptdata,2) % copy metadata headers
      SubjInfoAvg{1,i+ind-1} = ptdata{1,i};
    end

    for i=2:size(ptdata,1) % loop thru each subject in metadata file
      for j=2:size(ptdata,2) % loop thru each data column in metadata file
        SubjInfoAvg{i,j+ind-1} = ptdata{i,j};
      end
    end

    for i=1:size(SubjInfoAvg,1) % replace all empty cells to .
      for j=1:size(SubjInfoAvg,2)
        if isempty(SubjInfoAvg{i,j}), SubjInfoAvg{i,j} = '.'; end;
      end
    end

    write_cellarray_csv(outmetadatafnameavg,SubjInfoAvg);
    fprintf('%s: Complete files finished writing for project %s, normflag=%d\n',...
      mfilename,ProjID,normflag);
  end % end norm loop
end % end project loop
