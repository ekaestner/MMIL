function aseg_stats = fs_read_aseg_stats(subj,subjdir,roigroups);
%function aseg_stats = fs_read_aseg_stats(subj,[subjdir],[roigroups]);
%
% Required input:
%  subj: string specifying the subject name
%
% Optional input:
%  subjdir: subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  aseg_roigroups: struct array containing the following fields:
%    roiname: name of new ROI
%    roicode: new ROI code number
%    roicodes: vector of aseg ROI code numbers
%    {default: aseg_roigroups = fs_define_aseg_roigroups}
%     (includes 'WholeBrain', 'LatVentricles', and 'AllVentricles')
%
% Output:
%   aseg_stats is a struct array containing:
%     roiname
%     roicode
%     volume
%
% Created:  04/01/09 by Don Hagler
% Last Mod: 11/06/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aseg_stats = [];
if ~mmil_check_nargs(nargin,1), return; end;

if ~exist('subjdir','var') | isempty(subjdir)
  subjdir = getenv('SUBJECTS_DIR');
  if isempty(subjdir)
    error('SUBJECTS_DIR not defined as environment variable');
  end;
end;

aseg_fname = sprintf('%s/%s/stats/aseg.stats',subjdir,subj);
if ~exist(aseg_fname,'file')
  fprintf('%s: WARNING: aseg stats file %s not found\n',mfilename,aseg_fname);
  return;
end;

if ~exist('roigroups','var') | isempty(roigroups)
  roigroups = fs_define_aseg_roigroups;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read aseg file
if ~isempty(aseg_fname)
  try
    patt = '# cvs_version $Id: mri_segstats.c,v ';
    segstats_ver = 0;
    fid = fopen(aseg_fname);
    % Now get mri_segstats version that created aseg_fname.
    while 1,
       tline = fgetl(fid);
       if ~ischar(tline), break, end
       if ~isempty(strfind(tline,patt)),
          ver_str = strtok(tline(length(patt):end),' ');
          segstats_ver = str2double(regexp(ver_str,'^\d+\.\d+','match'));
          break;
       end
    end
    frewind(fid);

    tmp_stats = textscan(fid,'%d %d %d %f %s %f %f %f %f %f\n',...
      'commentstyle','#');
    for i=1:length(tmp_stats{1})
      aseg_stats(i).roiname = char(tmp_stats{5}{i});
      aseg_stats(i).roicode = double(tmp_stats{2}(i));
      aseg_stats(i).volume  = double(tmp_stats{4}(i));
    end;

    % add extra ROIs from stats file  
    ii=length(tmp_stats{1});

    if segstats_ver>=1.11 & segstats_ver<=1.33,
       asegROI.hdr_str = {'Brain Mask Volume', ...
                       'Brain Segmentation Volume', ...
                       'Intracranial Volume'};
       asegROI.name = {'BrainMaskVolume', 'BrainSegmentationVolume', 'IntracranialVolume'};
       asegROI.code = [20001 20002 20003];
    elseif segstats_ver>=1.75,
       asegROI.hdr_str = {'Left hemisphere cortical white matter volume', ...
                          'Left hemisphere cortical gray matter volume', ...
                          'Right hemisphere cortical white matter volume', ...
                          'Right hemisphere cortical gray matter volume', ...
                          'Intracranial Volume', 'Supratentorial volume', ...
                          'Subcortical gray matter volume'};
       asegROI.name = {'Left-Cortical-White-Matter','Left-Cortical-Gray-Matter', ...
                       'Right-Cortical-White-Matter','Right-Cortical-Gray-Matter', ...
                       'IntracranialVolume','SupratentorialVolume','SubCorticalGrayVolume'};
       asegROI.code = [2 3 41 42 20003 20009 20010];
    else
       error('mri_segstats version incompatible.');
    end

    frewind(fid); % go back to beginning of file
    while 1
       tmp = fgetl(fid);
       if ~ischar(tmp) || tmp(1)~='#', break; end

       for jj=1:length(asegROI.hdr_str),
         ind=strfind(tmp,sprintf(' %s,',asegROI.hdr_str{jj}));
         if length(ind)==1,
           ii=ii+1;
           aseg_stats(ii).roiname = asegROI.name{jj};
           aseg_stats(ii).roicode = asegROI.code(jj);
           [T,R]=strtok(tmp(ind:end),','); 
           [T,R]=strtok(R,',');
           aseg_stats(ii).volume = str2double(T); 
           if isnan(aseg_stats(ii).volume),
             fprintf('%s: WARNING: unable to read %s\n',mfilename,asegROI.name{jj});
           end
           break;
         end
       end
    end
    fclose(fid);

    % add extra ROIs from roigroups
    all_codes =  cell2mat({aseg_stats.roicode});
    for i=1:length(roigroups)
      [roicodes,ind_group,ind_all] = intersect(roigroups(i).roicodes,all_codes);
      if isempty(ind_all)
        error('invalid roicodes in roigroups');
      else
        volume = 0;
        for j=ind_all
          volume = volume + aseg_stats(j).volume;
        end;
        aseg_stats(end+1).roiname = roigroups(i).roiname;
        aseg_stats(end).roicode = roigroups(i).roicode;
        aseg_stats(end).roicodes = roicodes;
        aseg_stats(end).volume = volume;
      end;
    end
  catch
    fprintf('%s: ERROR: failed to read aseg stats file\n',mfilename);
  end;
end;
