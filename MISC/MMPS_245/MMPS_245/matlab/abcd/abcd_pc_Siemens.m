function abcd_pc_Siemens(metadata,n_files,outdir,csadata)
%function abcd_pc_Siemens(metadata,n_files,outdir)
%
% Created:  08/18/16 by Jose Teruel
% Prev Mod: 12/14/16 by Jose Teruel
% Last Mod: 06/22/17 by Feng Xue added coiltype/pixelbandwidth/channel
%

% NOTE: based on siemens_compliance_check, last mod 08/18/16, created by
% Jose Teruel

if ~mmil_check_nargs(nargin,3), return; end;
global version;
version = '0.0.4';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

expectedSequences = cell(5);
expectedSequences{1}= 'tfl3d1_16ns';
expectedSequences{2}= 'spc_200ns';
expectedSequences{3}= 'ep_b';
expectedSequences{4}= 'epse2d1';
expectedSequences{5}= 'epfSM2d1';
expectedSequences{6}= 'epfid2d1_90';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imageType = mmil_getfield(metadata,'ImageType');
sequenceType = mmil_getfield(metadata,'SequenceName');

if any(strfind(sequenceType,expectedSequences{1}))
  check_T1_Siemens(metadata,n_files,outdir);
elseif any(strfind(sequenceType,expectedSequences{2}))
  check_T2_Siemens(metadata,n_files,outdir);
elseif any(strfind(sequenceType,expectedSequences{3})) && ~any(strfind(imageType,'MOSAIC')) && ~any(strfind(imageType,'TRACE'))
  check_dMRI_FM_Siemens(metadata,n_files,outdir,csadata);
elseif any(strfind(sequenceType,expectedSequences{3})) && any(strfind(imageType,'MOSAIC'))
  check_dMRI_Siemens(metadata,n_files,outdir);
elseif any(strfind(sequenceType,expectedSequences{4}))
  check_fMRI_FM_Siemens(metadata,n_files,outdir,csadata);
elseif any(strfind(sequenceType,expectedSequences{5})) || any(strfind(sequenceType,expectedSequences{6}))
  check_fMRI_Siemens(metadata,n_files,outdir);    
else
  check_undefined(metadata,n_files,outdir);
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% T1
function check_T1_Siemens(metadata,n_files,outdir)
  seriesType = 'T1';
  ABCD_compliant = 'No';
  completed_series = 0;
  informationString ='SIEMENS T1-weighted; ';
  f_rows = 0; % flag rows
  f_pe_steps = 0; % flag phase encoding steps
  f_norm = 0; % flag NORM
  if ismember(mmil_getfield(metadata,'Rows'),[256])
    f_rows = 1;
    informationString = [informationString,...
      'Valid number of rows; '];
  else
    informationString = [informationString,...
      'Invalid number of rows; '];
  end
  if ismember(mmil_getfield(metadata,'NumberOfPhaseEncodingSteps'),[255])
    f_pe_steps = 1;
    informationString = [informationString,...
      'Valid number of phase encoding steps; '];
  else
    informationString = [informationString,...
      'Invalid number of phase encoding steps; '];
  end
  if any(strfind(mmil_getfield(metadata,'ImageType'), 'NORM'))  
    f_norm = 1;
    informationString = [informationString,...
      'Normalized T1 image; '];
    seriesType ='T1_NORM';
    ABCD_compliant = 'NA';
  elseif any(strfind(mmil_getfield(metadata,'SeriesDescription'), 'Collection')) 
    f_norm = 1;
    informationString = [informationString,...
      'MPR Collection image; '];
    seriesType ='T1_Collection';
    ABCD_compliant = 'NA';      
  else
    informationString = [informationString,...
      'Not normalized T1 image; '];
  end
   
  if (f_rows && f_pe_steps && ~f_norm)
    %seriesType = 'ABCD_T1';
    ABCD_compliant = 'Yes';
  end
  
  if n_files == 176
    completed_series = 1;
  end
  
  
  
  write_output_json(metadata,n_files,seriesType,informationString,ABCD_compliant,completed_series,outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% T2
function check_T2_Siemens(metadata,n_files,outdir)
  seriesType = 'T2';
  ABCD_compliant = 'No';
  informationString ='SIEMENS T2-weighted; '; 
  f_rows = 0; % flag Rows
  f_pe_steps = 0; % flag phase encoding steps
  f_norm = 0; % flag NORM
  completed_series = 0;
  if ismember(mmil_getfield(metadata,'Rows'),[256])
    f_rows = 1;
    informationString = [informationString,...
      'Valid number of rows; '];
  else
    informationString = [informationString,...
      'Invalid number of rows; '];
  end
  if ismember(mmil_getfield(metadata,'NumberOfPhaseEncodingSteps'),[227])
    f_pe_steps = 1;
    informationString = [informationString,...
      'Valid number of phase encoding steps; '];
  else
    informationString = [informationString,...
      'Invalid number of phase encoding steps; '];
  end
  if any(strfind(mmil_getfield(metadata,'ImageType'), 'NORM')) 
    f_norm = 1;
    informationString = [informationString,...
      'Normalized T2 image; '];
    seriesType ='T2_NORM';
    ABCD_compliant = 'NA';
  else
    informationString = [informationString,...
      'Not normalized T2 image; '];
  end

  if (f_rows && f_pe_steps && ~f_norm)
    %seriesType = 'ABCD_T2';
    ABCD_compliant = 'Yes';
  end
  
  if n_files == 176
    completed_series = 1;
  end
  
  write_output_json(metadata,n_files,seriesType,informationString,ABCD_compliant,completed_series,outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dMRI field map
function check_dMRI_FM_Siemens(metadata,n_files,outdir,csadata)
  seriesType = 'dMRI_FM';
  ABCD_compliant = 'No';
  informationString = 'SIEMENS Diffusion FM; ';
  f_rows = 0; % flag rows
  f_pe_steps = 0; % flag phase encoding steps
  completed_series = 0;
  pe_dir_pos = -1;
  %==================csa header====================%
  names = {csadata{1}.CSAImageHeaderInfo.name};
  idx = find(strcmp(names,'PhaseEncodingDirectionPositive'));
  if isempty(idx)
    fprintf('WARNING: PhaseEncodingDirectionPositive in CSA header not found\n');
  else
    pe_dir_pos = str2double(csadata{1}.CSAImageHeaderInfo(idx).item(1).val);
  end;
  %================================================%
   
  if pe_dir_pos==0
      seriesType = 'dMRI_FM_PA';
  elseif pe_dir_pos==1
      seriesType = 'dMRI_FM_AP';
  else
      informationString = [informationString,...
        'Phase encoding direction could not be established ;'];
  end
  if ismember(mmil_getfield(metadata,'Rows'),[140]) 
    f_rows = 1;
    informationString = [informationString,...
      'Valid number of rows; '];
  else
    informationString = [informationString,...
      'Invalid number of rows; '];
  end
  if ismember(mmil_getfield(metadata,'NumberOfPhaseEncodingSteps'),[105])  
    f_pe_steps = 1;
    informationString = [informationString,...
      'Valid number of phase encoding steps; '];
  else
    informationString = [informationString,...
      'Invalid number of phase encoding steps; '];
  end

  if (f_rows && f_pe_steps)
     ABCD_compliant = 'Yes';
  end
  
  if ismember(n_files, [81,82])
    completed_series = 1;
  end
      
  write_output_json(metadata,n_files,seriesType,informationString,ABCD_compliant,completed_series,outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dMRI
function check_dMRI_Siemens(metadata,n_files,outdir)
  seriesType = 'dMRI';
  ABCD_compliant = 'No';
  informationString = 'SIEMENS Diffusion; ';
  f_rows = 0; % flag Rows
  f_pe_steps = 0; % flag Phase encoding steps
  completed_series = 0;
  if ismember(mmil_getfield(metadata,'Rows'),[1260]) 
    f_rows = 1;
    informationString = [informationString,...
      'Valid number of rows; '];
  else
    informationString = [informationString,...
      'Invalid number of rows; '];
  end
  if ismember(mmil_getfield(metadata,'NumberOfPhaseEncodingSteps'),[105])  
    f_pe_steps = 1;
    informationString = [informationString,...
      'Valid number of phase encoding steps; '];
  else
    informationString = [informationString,...
      'Invalid number of phase encoding steps; '];
  end

  if (f_rows && f_pe_steps)
    %seriesType = 'ABCD_dMRI';
    ABCD_compliant = 'Yes';
  end
  
  if n_files == 103
     completed_series = 1;
  end  
  
  write_output_json(metadata,n_files,seriesType,informationString,ABCD_compliant,completed_series,outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fMRI field map
function check_fMRI_FM_Siemens(metadata,n_files,outdir,csadata)
  seriesType = 'fMRI_FM';
  ABCD_compliant = 'No';
  informationString = 'SIEMENS fMRI FM; ';
  % checks to perform
  f_rows = 0; % flag rows
  f_pe_steps = 0; % flag phase encoding steps
  f_norm = 0; % flag NORM
  completed_series = 0;
  pe_dir_pos = -1;
  %==================csa header====================%
  names = {csadata{1}.CSAImageHeaderInfo.name};
  idx = find(strcmp(names,'PhaseEncodingDirectionPositive'));
  if isempty(idx)
    fprintf('WARNING: PhaseEncodingDirectionPositive in CSA header not found\n');
  else
    pe_dir_pos = str2double(csadata{1}.CSAImageHeaderInfo(idx).item(1).val);
  end;
  %================================================%
  
  if pe_dir_pos==0
      seriesType = 'fMRI_FM_PA';
  elseif pe_dir_pos==1
      seriesType = 'fMRI_FM_AP';
  else
      informationString = [informationString,...
        'Phase encoding direction could not be established ;'];
  end
  
  if ismember(mmil_getfield(metadata,'Rows'),[90]) 
    f_rows = 1;
    informationString = [informationString,...
      'Valid number of rows; '];
  else
    informationString = [informationString,...
      'Invalid number of rows; '];
  end
  if ismember(mmil_getfield(metadata,'NumberOfPhaseEncodingSteps'),[90]) 
    f_pe_steps = 1;
    informationString = [informationString,...
      'Valid number of phase encoding steps; '];
  else
    informationString = [informationString,...
      'Invalid number of phase encoding steps; '];
  end
  if any(strfind(mmil_getfield(metadata,'ImageType'), 'NORM')) 
    f_norm = 1;
    informationString = [informationString,...
      'Normalized fMRI_FM image; '];
  else
    informationString = [informationString,...
      'Not normalized fMRI_FM image; '];
  end

  if (f_rows && f_pe_steps)
    ABCD_compliant = 'Yes';  
  end
  
  if n_files == 60
    completed_series = 1;
  end
  
  write_output_json(metadata,n_files,seriesType,informationString,ABCD_compliant,completed_series,outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fMRI
function check_fMRI_Siemens(metadata,n_files,outdir)
  seriesType = 'Undefined_fMRI';
  informationString = 'SIEMENS fMRI ';
  ABCD_compliant = 'No';
  f_rows = 0; % flag rows
  f_pe_steps = 0; % flag phase encoding steps
  f_time = 0; % flag for acquisition time
  completed_series = 0;
  
  if ismember(mmil_getfield(metadata,'Rows'),[720]) 
    f_rows = 1;
    informationString = [informationString,...
      'Valid number of rows; '];
  else
    informationString = [informationString,...
      'Invalid number of rows; '];
  end
  
  if ismember(mmil_getfield(metadata,'NumberOfPhaseEncodingSteps'),[90]) 
    f_pe_steps = 1;
    informationString = [informationString,...
      'Valid number of phase encoding steps; '];
  else
    informationString = [informationString,...
      'Invalid number of phase encoding steps; '];
  end
  
  % todo: check that this tag (Private_0051_100a) exists (not all have it)
  aquisition_dur = mmil_getfield(metadata,'Private_0051_100a');
  if any(strfind(aquisition_dur,'05:11'))
    ftime = 1;
    seriesType = 'rsfMRI';
    informationString = [informationString,...
      'Valid acquisition time; '];
    if n_files == 383
        completed_series = 1;
    end
    if (f_pe_steps && f_rows)
      ABCD_compliant = 'Yes';
    end
  elseif any(strfind(aquisition_dur,'05:22')) || any(strfind(aquisition_dur,'05:33'))
    ftime = 1;
    seriesType = 'fMRI_MID_task';
    informationString = [informationString,...
      'Valid acquisition time; '];
    if n_files == 411
        completed_series = 1;
    end
    if (f_pe_steps && f_rows)

      ABCD_compliant = 'Yes';
    end  
  elseif any(strfind(aquisition_dur,'06:00'))
    ftime = 1;
    seriesType = 'fMRI_SST_task';
    informationString = [informationString,...
      'Valid acquisition time; '];
    if n_files == 445
        completed_series = 1;
    end
    if (f_pe_steps && f_rows)
      ABCD_compliant = 'Yes';
    end
  elseif any(strfind(aquisition_dur,'05:00'))
    ftime = 1;
    seriesType = 'fMRI_nBack_task';
    informationString = [informationString,...
      'Valid acquisition time; '];
    if n_files == 370
        completed_series = 1;
    end
    if (f_pe_steps && f_rows)
      %seriesType = 'ABCD_rsfMRI';
      ABCD_compliant = 'Yes';
    end
  else
    informationString = [informationString,...
      'Invalid acquisition time'];
  end
  
  %%%If imaging parameters cannot distinguish the fMRI sequence we
  %%%look at Series Description. They will be not compliant but completness
  %%%is checked
  
  if strcmp(seriesType,'Undefined_fMRI')
    if any(strfind(mmil_getfield(metadata,'SeriesDescription'),'rsfMRI'))
        seriesType = 'rsfMRI';
        if n_files == 383
            completed_series = 1;
        end   
    elseif any(strfind(mmil_getfield(metadata,'SeriesDescription'),'Monetary'))
        seriesType = 'fMRI_MID_task';
        if n_files == 411
            completed_series = 1;
        end
    elseif any(strfind(mmil_getfield(metadata,'SeriesDescription'),'Emotional'))
        seriesType = 'fMRI_nBack_task';
        if n_files == 370
            completed_series = 1;
        end             
    elseif any(strfind(mmil_getfield(metadata,'SeriesDescription'),'Stop'))
        seriesType = 'fMRI_SST_task'; 
        if n_files == 445
            completed_series = 1;
        end            
    end
  end
  write_output_json(metadata,n_files,seriesType,informationString,ABCD_compliant,completed_series,outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% undefined MRI
function check_undefined(metadata,n_files,se_outdir)
  if regexp(lower(mmil_getfield(metadata,'SeriesDescription')),'loc')
    seriesType = 'Localizer';
    informationString = 'Series Localizer';
    ABCD_compliant = 'NA';
  elseif regexp(lower(mmil_getfield(metadata,'SeriesDescription')),'calibration')
    seriesType = 'Calibration';
    informationString = 'Calibration Series';
    ABCD_compliant = 'NA';
  elseif regexp(lower(mmil_getfield(metadata,'SeriesDescription')),'setter')
    seriesType = 'vNav_setter';
    informationString = 'vNav_setter';
    ABCD_compliant = 'NA';
  elseif regexp(lower(mmil_getfield(metadata,'SeriesDescription')),'pmu')
    seriesType = 'PMU';
    informationString = 'Physio';
    ABCD_compliant = 'NA';    
  else
    seriesType = 'Undefined';
    informationString = 'Intended ABCD series was not discovered';
    ABCD_compliant = 'No';   
  end
  completed_series = 0; %Default as it is irrelevant
  write_output_json(metadata,n_files,seriesType,informationString,ABCD_compliant,completed_series,se_outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% write output file
function write_output_json(meta,n_files,seriesType,informationString,ABCD_compliant,completed_series,outdir)
  global version;
    
  sName = strrep(mmil_getfield(meta,'SequenceName'), '/', '_');
  sName = strrep(sName, '\', '_');
  pID = strrep(mmil_getfield(meta,'PatientID'), '/', '_');
  iType = strrep(mmil_getfield(meta,'ImageType'), '/', '_');
  iType = strrep(iType, '\', '_');
  [Channel,CoilType] = abcd_get_coiltype(meta,'Siemens');
  
  info_PC = struct('SeriesDescription',mmil_getfield(meta,'SeriesDescription'),...
    'SeriesType',seriesType,...
    'ABCD_Compliant',ABCD_compliant,...
    'Completed', completed_series,...
    'AdditionalInfo',informationString,...
    'NumberOfFiles',n_files,...
    'ImagesInAcquisition','NA',...
    'AcquisitionTime',mmil_getfield(meta,'Private_0051_100a'),...
    'NumberOfTemporalPositions','NA',...
    'AcquisitionMatrix', mmil_getfield(meta,'AcquisitionMatrix'),...
    'Rows',mmil_getfield(meta,'Rows'),...
    'PercentPhaseFieldOfView','NA',...
    'NumberOfPhaseEncodingSteps',mmil_getfield(meta,'NumberOfPhaseEncodingSteps'),...
    'RepetitionTime',mmil_getfield(meta,'RepetitionTime'),...
    'EchoTime',mmil_getfield(meta,'EchoTime'),...
    'SeriesNumber',mmil_getfield(meta,'SeriesNumber'),...
    'Manufacturer',mmil_getfield(meta,'Manufacturer'),...
    'SequenceName',sName,...
    'PixelBandwidth',mmil_getfield(meta,'PixelBandwidth'),...
    'Channel',Channel,...
    'CoilType',CoilType,...
    'ImageType',iType,...
    'PatientID',pID,...
    'PatientFamilyName',mmil_getfield(meta.PatientName,'FamilyName'),...
    'StudyInstanceUID',mmil_getfield(meta,'StudyInstanceUID'),...
    'SeriesInstanceUID',mmil_getfield(meta,'SeriesInstanceUID'),...
    'StudyDate',mmil_getfield(meta,'StudyDate'),...
    'StudyTime',mmil_getfield(meta,'StudyTime'),...
    'SeriesTime',mmil_getfield(meta,'SeriesTime'),...
    'version',version);
  
  jsonfile = savejson('ProtocolCompliance', info_PC ,struct('FloatFormat','%.2f'));
  resultFile = 'QC_PC.json'; %% todo: outstem parameter
  file = fullfile(outdir, resultFile);
  jfileID = fopen(file,'w');
  fprintf(jfileID, jsonfile);
  fclose(jfileID);
return;


