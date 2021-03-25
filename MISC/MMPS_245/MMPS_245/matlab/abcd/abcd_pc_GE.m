function abcd_pc_GE(metadata,n_files,outdir)
%function abcd_pc_GE(metadata,n_files,outdir)
%
% Created:  08/18/16 by Jose Teruel
% Prev Mod: 06/22/17 by Feng Xue
% Last Mod: 09/13/17 by Don Hagler
%

% NOTE: based on ge_compliance_check, last mod 08/18/16, created by Jose Teruel

if ~mmil_check_nargs(nargin,3), return; end;
global version;
version = '0.0.5';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

expectedSequences = cell(5);
expectedSequences{1}= 'research/mprage_promo';
expectedSequences{2}= 'research/3dfse_promo';
expectedSequences{3}= 'research/ABCD/muxepi2';
expectedSequences{4}= 'research/ABCD/muxepi';
expectedSequences{5}= 'research/ABCD/epi_pepolar';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aMatrix = nonzeros(mmil_getfield(metadata,'AcquisitionMatrix'));
if isempty(aMatrix)
   check_undefined(metadata,n_files,se_outdir);
end
sequenceType = mmil_getfield(metadata,'Private_0019_109c');
if any(strfind(sequenceType,expectedSequences{1}))
  check_T1_GE(aMatrix,metadata,n_files,outdir);
elseif any(strfind(sequenceType,expectedSequences{2}))
  check_T2_GE(aMatrix,metadata,n_files,outdir);
elseif any(strfind(sequenceType,expectedSequences{3}))
  check_dMRI_GE(aMatrix,metadata,n_files,outdir);
elseif any(strfind(sequenceType,expectedSequences{4}))
  check_fMRI_GE(aMatrix,metadata,n_files,outdir);
elseif any(strfind(sequenceType,expectedSequences{5}))
  check_fMRI_FM_GE(aMatrix,metadata,n_files,outdir);
else
  check_undefined(metadata,n_files,outdir);
end


return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% T1
function check_T1_GE(aMatrix,metadata,n_files,outdir)
  seriesType = 'T1';
  ABCD_compliant = 'No';
  informationString ='GE T1-weighted; '; 
  cmatrix = 0; % resolution
  completed_series = 0;
  if ismember(aMatrix,[256])
    cmatrix = 1; %Correct resolution
    informationString = [informationString,...
      'Valid matrix size; '];
  else
    informationString = [informationString,...
      'Invalid matrix size; '];
  end
  
  if cmatrix
    %seriesType = 'ABCD_T1';
    ABCD_compliant = 'Yes';
  end
  
  %if (mmil_getfield(metadata,'ImagesInAcquisition') <= n_files)
  if (n_files == 208)
    completed_series = 1;
  end 
  
  write_output_json(metadata,n_files,seriesType,informationString,ABCD_compliant,completed_series,outdir);

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% T2
function check_T2_GE(aMatrix,metadata,n_files,outdir)
  seriesType = 'T2';
  ABCD_compliant = 'No';
  informationString ='GE T2-weighted; '; 
  cmatrix = 0; % resolution
  completed_series = 0;
  
  if ismember(aMatrix,[256])
    cmatrix = 1; % correct resolution
    informationString = [informationString,...
      'Valid matrix size; '];
  else
    informationString = [informationString,...
      'Invalid matrix size; '];
  end

  if (cmatrix)
    %seriesType = 'ABCD_T2';
    ABCD_compliant = 'Yes';
  end
  
  %if (mmil_getfield(metadata,'ImagesInAcquisition') <= n_files)
  if ismember(n_files, [208,204])
    completed_series = 1;
  end 
  
  
  write_output_json(metadata,n_files,seriesType,informationString,ABCD_compliant,completed_series,outdir);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fMRI field map
function check_fMRI_FM_GE(aMatrix,metadata,n_files,outdir)
  seriesType = 'fMRI_FM';
  ABCD_compliant = 'No';
  informationString ='GE fMRI_FM; '; 
  cmatrix = 0; % resolution
  cfiles = 0; % number of files
  completed_series = 0;
  nvolumes = 0; % number of volumes
  if ismember(aMatrix,[90])
    cmatrix = 1; % valid resolution
    informationString = [informationString,...
      'Valid matrix size; '];
  else
    informationString = [informationString,...
      'Invalid matrix size; '];
  end
  if ismember(mmil_getfield(metadata,'ImagesInAcquisition'),[120,180,240])
    cfiles = 1;
    informationString = [informationString,...
      'Valid number of images in acquisition; '];
  else
    informationString = [informationString,...
      'Invalid number of images in acquisition; '];
  end
  if (mmil_getfield(metadata,'Private_0019_10b3') == 2)
    nvolumes = 1;
    informationString = [informationString,...
      'Valid number of volumes; '];
  else
    informationString = [informationString,...
      'Invalid number of volumes; '];
  end
  if (cmatrix && cfiles && nvolumes)
    %seriesType = 'ABCD_fMRI_FM';
    ABCD_compliant = 'Yes';
  end
    
  %if (mmil_getfield(metadata,'ImagesInAcquisition') <= n_files)
  if ismember(n_files, [120,180,240])
    completed_series = 1;
  end 
  
  write_output_json(metadata,n_files,seriesType,informationString,ABCD_compliant,completed_series,outdir);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% diffusion field map or dMRI
function check_dMRI_GE(aMatrix,metadata,n_files,outdir)
  seriesType = 'Undefined_dMRI';
  ABCD_compliant = 'No';
  informationString = 'GE Diffusion; ';
  cfiles = 0;
  cmatrix = 0;
  completed_series = 0;
  if ismember(mmil_getfield(metadata,'ImagesInAcquisition'),[162,648])
    seriesType = 'dMRI_FM';
    cfiles = 1;
    informationString = [informationString,...
      'Field map;  Valid number of images in acquisition; '];
    if ismember(aMatrix,[140])
      cmatrix = 1;
      informationString = [informationString,...
        'Valid matrix size; '];
    else
      informationString = [informationString,...
        'Invalid matrix size; '];
    end
    if (cfiles && cmatrix)
      %seriesType = 'ABCD_Diffusion_FM';
      ABCD_compliant = 'Yes';
    end
    
    %if (mmil_getfield(metadata,'ImagesInAcquisition') <= n_files)
    if ismember(n_files, [162,648])
      completed_series = 1;
    end
    
  elseif ismember(mmil_getfield(metadata,'ImagesInAcquisition'),[8343,8424]) % Number of files after recon (2916 if not)
    seriesType = 'dMRI';
    cfiles = 1;
    informationString = [informationString,...
      'DTI;  Valid number of images in acquisition; '];
    if ismember(aMatrix,[140])
      cmatrix = 1;
      informationString = [informationString,...
        'Valid matrix size; '];
    else
      informationString = [informationString,...
        'Invalid matrix size; '];
    end

    if (cfiles && cmatrix)
      %seriesType = 'ABCD_DTI';
      ABCD_compliant = 'Yes';
    end
    
    %if (mmil_getfield(metadata,'ImagesInAcquisition') <= n_files)
    if ismember(n_files, [8343,8424])
      completed_series = 1;
    end
     
  else
    cfiles = 0;
    informationString = [informationString,...
      'Invalid number of images in acquisition; '];
  end  
  
  %%%If imaging parameters cannot distinguish between dMRI and dMRI_FM we
  %%%look at Series Description. They will be not compliant but completness
  %%%is checked
    
  if strcmp(seriesType,'Undefined_dMRI')
      if any(strfind(mmil_getfield(metadata,'SeriesDescription'),'Map'))
        seriesType = 'dMRI_FM';
        %if (mmil_getfield(metadata,'ImagesInAcquisition') <= n_files)
        if ismember(n_files, [162,648])
          completed_series = 1;
        end      
      elseif any(strfind(mmil_getfield(metadata,'SeriesDescription'),'DTI'))
        seriesType = 'dMRI';
        %if (mmil_getfield(metadata,'ImagesInAcquisition') <= n_files)
        if ismember(n_files, [8343,8424])
          completed_series = 1;
        end        
      end
  end
    
  write_output_json(metadata,n_files,seriesType,informationString,ABCD_compliant,completed_series,outdir);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fMRI
function check_fMRI_GE(aMatrix,metadata,n_files,outdir)
  seriesType = 'Undefined_fMRI';
  informationString = 'GE fMRI ';
  ABCD_compliant = 'No';
  cfiles = 0;
  cmatrix = 0;
  completed_series = 0;
  aquisition_dur = mmil_getfield(metadata,'Private_0019_105a');
  if ismember(aquisition_dur,[312800384])
    seriesType = 'rsfMRI';
    cfiles = 1;
    informationString = [informationString,...
      'Resting State; Acquisition time within range; '];
    if ismember(aMatrix,[90])
      cmatrix = 1;
      informationString = [informationString,...
        'Valid matrix size; '];
    else
      informationString = [informationString,...
        'Invalid matrix size; '];
    end
    %if (mmil_getfield(metadata,'ImagesInAcquisition') <= n_files)
    if ismember(n_files, [22800,23460])
        completed_series = 1;
    end     
    if (cfiles && cmatrix)
      %seriesType = 'ABCD_rsfMRI';
      ABCD_compliant = 'Yes';
    end
  elseif ismember(aquisition_dur,[302400384])
    seriesType = 'fMRI_nBack_task';
    cfiles = 1;
    informationString = [informationString,...
      'fMRI nBack task; Acquisition time within range; '];
    if ismember(aMatrix,[90])
      cmatrix = 1;
      informationString = [informationString,...
        'Valid matrix size; '];
    else
      informationString = [informationString,...
        'Invalid matrix size; '];
    end
    %if (mmil_getfield(metadata,'ImagesInAcquisition') <= n_files)
    if ismember(n_files, [22020,22680])
        completed_series = 1;
    end 
    if (cfiles && cmatrix)
      %seriesType = 'ABCD_fMRI_nBack_task';
      ABCD_compliant = 'Yes';
    end 
  elseif ismember(aquisition_dur,[362400448])
    seriesType = 'fMRI_SST_task';
    cfiles = 1;
    informationString = [informationString,...
      'fMRI SST task; Acquisition time within range; '];
    if ismember(aMatrix,[90])
      cmatrix = 1;
      informationString = [informationString,...
        'Valid matrix size; '];
    else
      informationString = [informationString,...
        'Invalid matrix size; '];
    end
    %if (mmil_getfield(metadata,'ImagesInAcquisition') <= n_files)
    if ismember(n_files, [26520,27180])
        completed_series = 1;
    end 
    if (cfiles && cmatrix)
      %seriesType = 'ABCD_fMRI_SST_task';
      ABCD_compliant = 'Yes';
    end  
  elseif ismember(aquisition_dur,[335200416])
    seriesType = 'fMRI_MID_task';
    cfiles = 1;
    informationString = [informationString,...
      'fMRI MID task; Acquisition time within range; '];
    if ismember(aMatrix,[90])
      cmatrix = 1;
      informationString = [informationString,...
        'Valid matrix size; '];
    else
      informationString = [informationString,...
        'Invalid matrix size; '];
    end
    %if (mmil_getfield(metadata,'ImagesInAcquisition') <= n_files)
    if ismember(n_files, [24480,25140])
        completed_series = 1;
    end 
    if (cfiles && cmatrix)
      %seriesType = 'ABCD_fMRI_MID_task';
      ABCD_compliant = 'Yes';
    end  
  else
    informationString = [informationString,...
      'Invalid acquisition time; '];
  end
  
  %%%If imaging parameters cannot distinguish the fMRI sequence we
  %%%look at Series Description. They will be not compliant but completness
  %%%is checked
  
    if strcmp(seriesType,'Undefined_fMRI')
        if any(strfind(mmil_getfield(metadata,'SeriesDescription'),'rsfMRI'))
            seriesType = 'rsfMRI';
            %if (mmil_getfield(metadata,'ImagesInAcquisition') <= n_files)
            %if (mmil_getfield(metadata,'ImagesInAcquisition') <= n_files)
            if ismember(n_files, [22800,23460])
                completed_series = 1;
            end      
        elseif any(strfind(mmil_getfield(metadata,'SeriesDescription'),'Monetary'))
            seriesType = 'fMRI_MID_task';
            %if (mmil_getfield(metadata,'ImagesInAcquisition') <= n_files)
            if ismember(n_files, [24480,25140])
                completed_series = 1;
            end  
        elseif any(strfind(mmil_getfield(metadata,'SeriesDescription'),'Emotional'))
            seriesType = 'fMRI_nBack_task';
             %if (mmil_getfield(metadata,'ImagesInAcquisition') <= n_files)
             if ismember(n_files, [22020,22680])
                completed_series = 1;
            end              
        elseif any(strfind(mmil_getfield(metadata,'SeriesDescription'),'Stop'))
            seriesType = 'fMRI_SST_task'; 
            %if (mmil_getfield(metadata,'ImagesInAcquisition') <= n_files)
            if ismember(n_files, [26520,27180])
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

function write_output_json(meta,n_files,seriesType,informationString,ABCD_compliant,completed_series,outdir)
  global version;
  

  sName = strrep(mmil_getfield(meta,'Private_0019_109c'), '/', '_');
  sName = strrep(sName, '\', '_');
  pID = strrep(mmil_getfield(meta,'PatientID'), '/', '_');
  iType = strrep(mmil_getfield(meta,'ImageType'), '/', '_');
  iType = strrep(iType, '\', '_');  
  [Channel,CoilType] = abcd_get_coiltype(meta,'GE');


  info_PC = struct('SeriesDescription',mmil_getfield(meta,'SeriesDescription'),...
    'SeriesType',seriesType,...
    'ABCD_Compliant',ABCD_compliant,...
    'Completed', completed_series,...
    'AdditionalInfo',informationString,...
    'NumberOfFiles',n_files,...
    'ImagesInAcquisition',mmil_getfield(meta,'ImagesInAcquisition'),...
    'AcquisitionTime',mmil_getfield(meta,'Private_0019_105a'),...
    'NumberOfTemporalPositions',mmil_getfield(meta,'NumberOfTemporalPositions'),...
    'AcquisitionMatrix', mmil_getfield(meta,'AcquisitionMatrix'),...
    'Rows',mmil_getfield(meta,'Rows'),...
    'PercentPhaseFieldOfView',mmil_getfield(meta,'PercentPhaseFieldOfView'),...
    'NumberOfPhaseEncodingSteps','NA',...
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
  resultFile = 'QC_PC.json'; %% todo: outstem parameter
  file = fullfile(outdir, resultFile);
  jfileID = fopen(file,'w');
  fprintf(jfileID, jsonfile);
  fclose(jfileID);
 
return;

