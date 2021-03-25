function abcd_pc_Philips(metadata,n_files,se_outdir)
%function abcd_pc_Philips(metadata,n_files,se_outdir)
%
% Created:  08/30/16 by Don Hagler
% Prev Mod: 11/23/16 by Jose Teruel
% Last Mod: 06/22/17 by Feng Xue added coiltype/pixelbandwidth/channel
%

% NOTE: based on siemens_compliance_check, last mod 08/18/16, created by
% Jose Teruel

if ~mmil_check_nargs(nargin,3), return; end;
global version;
version = '0.0.4';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

expectedSequences = cell(5);
expectedSequences{1}= 'T1TFE';
expectedSequences{2}= 'TSE';
expectedSequences{3}= 'DwiSE';
expectedSequences{4}= 'SEEPI';
expectedSequences{5}= 'FEEPI';
expectedSequences{6}= 'TFE';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metadata.AcquisitionMatrix = nonzeros(mmil_getfield(metadata,'AcquisitionMatrix'));
if isempty(metadata.AcquisitionMatrix)
   check_undefined(metadata,n_files,se_outdir);
end
sequenceType = mmil_getfield(metadata,'Private_2001_1020');
if any(strfind(sequenceType,expectedSequences{1}))
  check_T1_Philips(metadata,n_files,se_outdir);
elseif any(strfind(sequenceType,expectedSequences{6}))
  check_T1_Philips(metadata,n_files,se_outdir);
elseif any(strfind(sequenceType,expectedSequences{2}))
  check_T2_Philips(metadata,n_files,se_outdir);
elseif any(strfind(sequenceType,expectedSequences{3}))
  check_dMRI_Philips(metadata,n_files,se_outdir);
elseif any(strfind(sequenceType,expectedSequences{4}))
  check_FM_Philips(metadata,n_files,se_outdir);
elseif any(strfind(sequenceType,expectedSequences{5}))
  check_fMRI_Philips(metadata,n_files,se_outdir);    
else
  check_undefined(metadata,n_files,se_outdir);
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% T1
function check_T1_Philips(metadata,n_files,se_outdir)
  seriesType = 'Undefined';
  ABCD_compliant = 'No';
  informationString ='PHILIPS T1-weighted; '; 
  f_scanningsequence = 0; % flag scanning sequence
  f_sequencevariant = 0; % flag sequence variant
  f_noprojection = 0; %flag no projection image
  f_matrix = 0; % flag matrix
  completed_series = 0;
  if ismember(unique(mmil_getfield(metadata,'AcquisitionMatrix')),[256])
    f_matrix = 1;
    informationString = [informationString,...
      'Valid acquisition matrix; '];
  else
    informationString = [informationString,...
      'Invalid acquisition matrix; '];
  end
  if strcmp(mmil_getfield(metadata,'ScanningSequence'),'GR')
    f_scanningsequence = 1;
    informationString = [informationString,...
      'Correct scanning sequence; '];
  else
    informationString = [informationString,...
      'Incorrect scanning sequence; '];
  end
  if strcmp(mmil_getfield(metadata,'SequenceVariant'),'MP')
    f_sequencevariant = 1;
    informationString = [informationString,...
      'Correct sequence variant; '];
  else
    informationString = [informationString,...
      'Incorrect sequence variant; '];
  end
  if any(strfind(mmil_getfield(metadata,'ImageType'), 'PROJECTION IMAGE'))  
    informationString = [informationString,...
      'Projection image; '];
    ABCD_compliant = 'NA';
  else
    f_noprojection = 1;
    informationString = [informationString,...
      'Not projection image; '];
  end

  if (f_matrix && f_sequencevariant && f_scanningsequence && f_noprojection)
    seriesType = 'T1';
    ABCD_compliant = 'Yes';
  end
  
  if n_files >= 225
    completed_series = 1;
  end
  
  %This is a particular case as T1 for Philips use the same exact sequence
  %for localizers and other non-relevant series. In case the series is
  %Undefined based on sequence parameters we check the SeriesDescription
  
  if strcmp(seriesType,'Undefined')
    if strcmp(mmil_getfield(metadata,'SeriesDescription'),'3D T1')
        seriesType = 'T1';
    elseif (any(strfind(mmil_getfield(metadata,'SeriesDescription'),'Loc')) ||...
      any(strfind(mmil_getfield(metadata,'SeriesDescription'),'loc')) ||...
      any(strfind(mmil_getfield(metadata,'SeriesDescription'),'LOC')))
      seriesType = 'Localizer';
      informationString = 'Series Localizer';
      ABCD_compliant = 'NA';
    elseif (any(strfind(mmil_getfield(metadata,'SeriesDescription'),'SMARTPLAN')))
      seriesType = 'SMARTPLAN alignment';
      informationString = 'Philips SMARTPLAN alignment';
      ABCD_compliant = 'NA'; 
    end
  end
  write_output_json(metadata,n_files,seriesType,informationString,ABCD_compliant,completed_series,se_outdir);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% T2
function check_T2_Philips(metadata,n_files,se_outdir)
  seriesType = 'T2';
  ABCD_compliant = 'No';
  informationString ='PHILIPS T2-weighted; '; 
  f_scanningsequence = 0; % flag scanning sequence
  f_sequencevariant = 0; % flag sequence variant
  f_noprojection = 0; %flag no projection image
  f_matrix = 0; % flag matrix
  completed_series = 0;
  
  if ismember(unique(mmil_getfield(metadata,'AcquisitionMatrix')),[256])
    f_matrix = 1;
    informationString = [informationString,...
      'Valid acquisition matrix; '];
  else
    informationString = [informationString,...
      'Invalid acquisition matrix; '];
  end
  if strcmp(mmil_getfield(metadata,'ScanningSequence'),'SE')
    f_scanningsequence = 1;
    informationString = [informationString,...
      'Correct scanning sequence; '];
  else
    informationString = [informationString,...
      'Incorrect scanning sequence; '];
  end
  if strcmp(mmil_getfield(metadata,'SequenceVariant'),'SK')
    f_sequencevariant = 1;
    informationString = [informationString,...
      'Correct sequence variant; '];
  else
    informationString = [informationString,...
      'Incorrect sequence variant; '];
  end
  if any(strfind(mmil_getfield(metadata,'ImageType'), 'PROJECTION IMAGE'))  
    informationString = [informationString,...
      'Projection image; '];
    seriesType = 'T2_PROJECTION';
  else
    f_noprojection = 1;
    informationString = [informationString,...
      'Not projection image; '];
    ABCD_compliant = 'NA';
  end

  if (f_matrix && f_sequencevariant && f_scanningsequence && f_noprojection)
    ABCD_compliant = 'Yes';
  end
  
  if n_files >= 256
    completed_series = 1;
  end  
  
  write_output_json(metadata,n_files,seriesType,informationString,ABCD_compliant,completed_series,se_outdir);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% diffusion MRI or diffusion FM
function check_dMRI_Philips(metadata,n_files,se_outdir)
  seriesType = 'Undefined_dMRI';
  ABCD_compliant = 'No';
  informationString = 'PHILIPS Diffusion; ';
  f_matrix = 0; % flag matrix
  completed_series = 0;
  
  if (metadata.AcquisitionMatrix(1) == 140 && metadata.AcquisitionMatrix(2) == 141)  
    f_matrix = 1;
    informationString = [informationString,...
      'Valid acquisition matrix; '];
  else
    informationString = [informationString,...
      'Invalid acquisition matrix; '];
  end

  if (f_matrix)    
    if any(strfind(mmil_getfield(metadata,'SeriesDescription'),'Fieldmap')) 
        if any(strfind(mmil_getfield(metadata,'SeriesDescription'),'P'))
            ABCD_compliant = 'Yes';
            seriesType = 'dMRI_FM_PA';
            if ismember(n_files,[80,81,240]) %Minimun required for a Philips site (80, 81, 82)
                completed_series = 1;
            end
        elseif any(strfind(mmil_getfield(metadata,'SeriesDescription'),'A'))
            ABCD_compliant = 'Yes';
            seriesType = 'dMRI_FM_AP';
            if ismember(n_files,[80,81,240])
                completed_series = 1;
            end
        else
            seriesType = 'dMRI_FM';
            informationString = [informationString,...
            'Phase encoding direction could not be established ;'];
        end        
    else
        seriesType = 'dMRI'; 
        ABCD_compliant = 'Yes';
        if ismember(n_files,[4080,4131]) % Either 4080 or 4131
            completed_series = 1;
        end
    end
  end
   
  write_output_json(metadata,n_files,seriesType,informationString,ABCD_compliant,completed_series,se_outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Field map (FM or dMRI)
function check_FM_Philips(metadata,n_files,se_outdir)
  seriesType = 'Undefined_FM';
  ABCD_compliant = 'No';
  informationString = 'PHILIPS Field Map; ';
  f_matrix = 0; % flag rows
  completed_series = 0;
  
  if (metadata.AcquisitionMatrix(1) == 140 && metadata.AcquisitionMatrix(2) == 141)  
    seriesType = 'dMRI FM';
    informationString = [informationString,'Diffusion; Valid acquisition matrix; '];
    f_matrix = 1;
    if any(strfind(mmil_getfield(metadata,'SeriesDescription'),'P'))
        seriesType = 'dMRI_FM_PA';
    elseif any(strfind(mmil_getfield(metadata,'SeriesDescription'),'A'))
        seriesType = 'dMRI_FM_AP';
    else
      informationString = [informationString,...
        'Phase encoding direction could not be established ;'];
    end    
    if ismember(n_files,[80,81,240])
        completed_series = 1;
    end 
    if (f_matrix)
        ABCD_compliant = 'Yes';
    end
  
  elseif (metadata.AcquisitionMatrix(1) == 92 && metadata.AcquisitionMatrix(2) == 89)
    seriesType = 'fMRI_FM';
    informationString = [informationString,'fMRI; Valid acquisition matrix; '];
    f_matrix = 1;
    if any(strfind(mmil_getfield(metadata,'SeriesDescription'),'P'))
        seriesType = 'fMRI_FM_PA';
    elseif any(strfind(mmil_getfield(metadata,'SeriesDescription'),'A'))
        seriesType = 'fMRI_FM_AP';
    else
      informationString = [informationString,...
        'Phase encoding direction could not be established ;'];        
    end
    if ismember(n_files,[60,180])
        completed_series = 1;
    end 
    if (f_matrix)
        ABCD_compliant = 'Yes';
    end  
  else
    informationString = [informationString,...
      'Invalid acquisition matrix'];
  end  
  
  %%% If the sequence parameters cannot differentiate the inteded field map
  %%% we test with the SeriesDescription. There will be non-compliant but
  %%% completeness is checked
  
  if strcmp(seriesType,'Undefined_FM')
      if any(strfind(mmil_getfield(metadata,'SeriesDescription'),'fMRI Fieldmap P'))
         seriesType = 'fMRI_FM_PA'; 
         if ismember(n_files,[60,180])
            completed_series = 1;
         end
      elseif any(strfind(mmil_getfield(metadata,'SeriesDescription'),'fMRI Fieldmap A'))
         seriesType = 'fMRI_FM_AP'; 
         if ismember(n_files,[60,180])
            completed_series = 1;
         end
      elseif any(strfind(mmil_getfield(metadata,'SeriesDescription'),'DTI Fieldmap A'))
         seriesType = 'dMRI_FM_AP'; 
         if ismember(n_files,[80,81,240])
            completed_series = 1;
         end
      elseif any(strfind(mmil_getfield(metadata,'SeriesDescription'),'DTI Fieldmap P'))
         seriesType = 'dMRI_FM_PA'; 
         if ismember(n_files,[80,81,240])
            completed_series = 1;
         end            
      end
  end
  
  write_output_json(metadata,n_files,seriesType,informationString,ABCD_compliant,completed_series,se_outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function check_fMRI_Philips(metadata,n_files,se_outdir)
  seriesType = 'fMRI';
  informationString = 'PHILIPS fMRI; ';
  ABCD_compliant = 'No';
  f_matrix = 0; % flag rows
  ftempPos = 0; % flag for number of temporal positions
  completed_series = 0;
    
  if (metadata.AcquisitionMatrix(1) == 92 && metadata.AcquisitionMatrix(2) == 89)  
    f_matrix = 1;
    informationString = [informationString,...
      'Valid acquisition matrix; '];
  else
    informationString = [informationString,...
      'Invalid acquisition matrix; '];
  end 
  ntp = mmil_getfield(metadata,'NumberOfTemporalPositions');
  if ismember(ntp,[375,383])
    ftempPos = 1;
    seriesType = 'rsfMRI';
    informationString = [informationString,'Valid number of temporal positions; '];
    
    if n_files >= 22980
        completed_series = 1;
    end
    
    if (f_matrix && ftempPos)
      %seriesType = 'ABCD_rsfMRI';
      ABCD_compliant = 'Yes';
    end
    
  elseif ismember(ntp,[437,445])
    ftempPos = 1;
    seriesType = 'fMRI_SST_task';
    informationString = [informationString,'Valid number of temporal positions; '];
    
    if n_files >= 26700
        completed_series = 1;
    end
    
    if (f_matrix && ftempPos)
      %seriesType = 'ABCD_fMRI_SST_task';
      ABCD_compliant = 'Yes';
    end
  elseif ismember(ntp,[350,370])
    ftempPos = 1;
    seriesType = 'fMRI_nBack_task';
    informationString = [informationString,'Valid number of temporal positions; '];
    
    if n_files >= 22200
        completed_series = 1;
    end
    
    if (f_matrix && ftempPos)
      %seriesType = 'ABCD_fMRI_nBack_task';
      ABCD_compliant = 'Yes';
    end
  
  elseif ismember(ntp,[403,411])
    ftempPos = 1;
    seriesType = 'fMRI_MID_task';
    informationString = [informationString,'Valid number of temporal positions; '];

    if n_files >= 24660
        completed_series = 1;
    end
    
    if (f_matrix && ftempPos)
      %seriesType = 'ABCD_fMRI_MID_task';
      ABCD_compliant = 'Yes';
    end   
  else
    informationString = [informationString,...
      'Invalid number of temporal positions'];
  end 
  
  %%%If imaging parameters cannot distinguish the fMRI sequence we
  %%%look at Series Description. They will be not compliant but completness
  %%%is checked
  
  if strcmp(seriesType,'Undefined_fMRI')
    if any(strfind(mmil_getfield(metadata,'SeriesDescription'),'rsfMRI'))
        seriesType = 'rsfMRI';
        if n_files >= 22980
            completed_series = 1;
        end   
    elseif any(strfind(mmil_getfield(metadata,'SeriesDescription'),'MID'))
        seriesType = 'fMRI_MID_task';
        if n_files >= 24660
            completed_series = 1;
        end
    elseif any(strfind(mmil_getfield(metadata,'SeriesDescription'),'Emotion'))
        seriesType = 'fMRI_nBack_task';
        if n_files >= 22200
            completed_series = 1;
        end             
    elseif any(strfind(mmil_getfield(metadata,'SeriesDescription'),'SST'))
        seriesType = 'fMRI_SST_task'; 
        if n_files >= 26700
            completed_series = 1;
        end            
    end
  end
 
  write_output_json(metadata,n_files,seriesType,informationString,ABCD_compliant,completed_series,se_outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% undefined MRI
function check_undefined(metadata,n_files,se_outdir)
  if regexp(lower(mmil_getfield(metadata,'SeriesDescription')),'loc')
    seriesType = 'Localizer';
    informationString = 'Series Localizer';
    ABCD_compliant = 'NA';
  elseif regexp(lower(mmil_getfield(metadata,'SeriesDescription')),'smartplan')
    seriesType = 'SMARTPLAN alignment';
    informationString = 'Philips SMARTPLAN alignment';
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

function write_output_json(meta,n_files,seriesType,informationString,ABCD_compliant,completed_series,se_outdir)
  global version;
 
  sName = strrep(mmil_getfield(meta,'Private_2001_1020'), '/', '_');
  sName = strrep(sName, '\', '_');
  pID = strrep(mmil_getfield(meta,'PatientID'), '/', '_');
  iType = strrep(mmil_getfield(meta,'ImageType'), '/', '_');
  iType = strrep(iType, '\', '_');
  [Channel,CoilType] = abcd_get_coiltype(meta,'Philips');
  
  
  info_PC = struct('SeriesDescription',mmil_getfield(meta,'SeriesDescription'),...
    'SeriesType',seriesType,...
    'ABCD_Compliant',ABCD_compliant,...
    'Completed', completed_series,...
    'AdditionalInfo',informationString,...
    'NumberOfFiles',n_files,...
    'ImagesInAcquisition','NA',...
    'AcquisitionTime','NA',...
    'NumberOfTemporalPositions',mmil_getfield(meta,'NumberOfTemporalPositions'),...
    'AcquisitionMatrix', mmil_getfield(meta,'AcquisitionMatrix'),...
    'Rows',mmil_getfield(meta,'Rows'),...
    'PercentPhaseFieldOfView','NA',...
    'NumberOfPhaseEncodingSteps',mmil_getfield(meta,'NumberOfPhaseEncodingSteps'),...
    'RepetitionTime',mmil_getfield(meta,'RepetitionTime'),...
    'EchoTime',mmil_getfield(meta,'EchoTime'),...
    'SeriesNumber',mmil_getfield(meta,'SeriesNumber'),...
    'Manufacturer',mmil_getfield(meta,'Manufacturer'),...
    'SequenceName',sName,...
    'ImageType',iType,...
    'PixelBandwidth',mmil_getfield(meta,'PixelBandwidth'),...
    'Channel',Channel,...
    'CoilType',CoilType,...
    'PatientID',pID,...
    'PatientFamilyName',mmil_getfield(meta.PatientName,'FamilyName'),...
    'StudyInstanceUID',mmil_getfield(meta,'StudyInstanceUID'),...
    'SeriesInstanceUID',mmil_getfield(meta,'SeriesInstanceUID'),...
    'StudyDate',mmil_getfield(meta,'StudyDate'),...
    'StudyTime',mmil_getfield(meta,'StudyTime'),...
    'SeriesTime',mmil_getfield(meta,'SeriesTime'),...
    'version',version);
  
  jsonfile = savejson('ProtocolCompliance', info_PC ,struct('FloatFormat','%.2f'));
  resultFile = 'QC_PC.json';
  file = fullfile(se_outdir, resultFile);
  jfileID = fopen(file,'w');
  fprintf(jfileID, jsonfile);
  fclose(jfileID);
return;
