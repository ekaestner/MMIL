function abcd_pqa_summarize(varargin)
%function abcd_pqa_summarize(varargin)
%abcd_pqa_summarize Summary of this function goes here
%   Detailed explanation goes here
%%%%%%%%%%% Variables
%
%
%  Created:                by Jose Teruel
%  Last Mod by: 09/28/2017 by Feng Xue
%%%%%%%%%%% Variables

% check input parameters
parms = check_input(varargin);

load(parms.fname_siteinfo);
%%%%%%%%%%% Requirements

%%%%%%%% Usage and target folder %%%%%%%%%%%%%%%

version = 'ani';

%%%%%%%%%%%%%%%% Get all folders %%%%%%%%%%%%%%%%%

inputFoldersList = dir(parms.indir);

cmd = sprintf('find %s -type f -name "*.json" > %s/jsonlist.txt',parms.indir,parms.outdir);
unix(cmd);

%%%%%%%%%%% Read file %%%%%%%%%%%%%%%

fid=fopen(fullfile(parms.outdir,'jsonlist.txt'));
tline = fgetl(fid);
tlines = cell(0,1);
while ischar(tline)
    tlines{end+1,1} = tline;
    tline = fgetl(fid);
end
fclose(fid);

%%%%%%%%% Open each json file %%%%%%%%

n_series = length(tlines);

index = 0;
counter = 1;

while index<n_series
    index = index+1;
    currentME = loadjson(tlines{index});
    [~,st_folder] = fileparts(fileparts(tlines{index}));
    [site, phantom] = lookupSP(st_folder, currentME.fBIRN_Phantom_QA.SeriesInfo.ScannerSerialNumber,pqa_sites);
    StudyInstaceUID{counter,1} = currentME.fBIRN_Phantom_QA.SeriesInfo.StudyInstanceUID;
    StudyDate{counter,1} = currentME.fBIRN_Phantom_QA.SeriesInfo.StudyDate;
    StudyTime{counter,1} = currentME.fBIRN_Phantom_QA.SeriesInfo.StudyTime;
    Site{counter,1} = site;
    PhantomCode{counter,1} = phantom;
    Manufacturer{counter,1} = currentME.fBIRN_Phantom_QA.SeriesInfo.Manufacturer;
    Coil{counter,1} = currentME.fBIRN_Phantom_QA.SeriesInfo.Coil;
    SeriesTime{counter,1} = currentME.fBIRN_Phantom_QA.SeriesInfo.SeriesTime;
    SeriesNumber{counter,1} = uint16(currentME.fBIRN_Phantom_QA.SeriesInfo.SeriesNumber);
    SeriesDescription{counter,1} = currentME.fBIRN_Phantom_QA.SeriesInfo.SeriesDescription;
    RepetitionTime{counter,1} = currentME.fBIRN_Phantom_QA.SeriesInfo.RepetitionTime;
    EchoTime{counter,1} = currentME.fBIRN_Phantom_QA.SeriesInfo.EchoTime;
    FlipAngle{counter,1} = currentME.fBIRN_Phantom_QA.SeriesInfo.FlipAngle;
    Mean{counter,1} = currentME.fBIRN_Phantom_QA.QA_metrics.mean;
    StandardDev{counter,1} = currentME.fBIRN_Phantom_QA.QA_metrics.sd;
    SNR{counter,1} = currentME.fBIRN_Phantom_QA.QA_metrics.snr;
    SFNR{counter,1} = currentME.fBIRN_Phantom_QA.QA_metrics.sfnr;
    RMS{counter,1} = currentME.fBIRN_Phantom_QA.QA_metrics.rms;
    TemporalDrift{counter,1} = currentME.fBIRN_Phantom_QA.QA_metrics.temp_drift;
    TemporalDrift_per_minute{counter,1} = currentME.fBIRN_Phantom_QA.QA_metrics.temp_drift_per_minute;
    Max_Absolute_Temporal_Drift{counter,1} = currentME.fBIRN_Phantom_QA.QA_metrics.max_temp_drift;
    RDC{counter,1} = currentME.fBIRN_Phantom_QA.QA_metrics.rdc;
    FWHM_x{counter,1} = currentME.fBIRN_Phantom_QA.QA_metrics.FWHM_x;
    FWHM_y{counter,1} = currentME.fBIRN_Phantom_QA.QA_metrics.FWHM_y;
    FWHM_z{counter,1} = currentME.fBIRN_Phantom_QA.QA_metrics.FWHM_z;
    MeanGhost{counter,1} = currentME.fBIRN_Phantom_QA.QA_metrics.mean_ghost;
    TopGhost{counter,1} = currentME.fBIRN_Phantom_QA.QA_metrics.top_ghost;
    SpatialDrift_PE{counter,1} = currentME.fBIRN_Phantom_QA.QA_metrics.spatial_drift_pe;
    ImageFrequency{counter,1} = currentME.fBIRN_Phantom_QA.QA_metrics.image_frequency;
    ScannerSerialNumber{counter,1} = currentME.fBIRN_Phantom_QA.SeriesInfo.ScannerSerialNumber;
    OperatingSystemVersion{counter,1} = currentME.fBIRN_Phantom_QA.SeriesInfo.OperatingSystemVersion;
    if isfield(currentME.fBIRN_Phantom_QA.SeriesInfo,'OSLevel');
        OSLevel{counter,1} = currentME.fBIRN_Phantom_QA.SeriesInfo.OSLevel;
    else
        OSLevel{counter,1} = '';
    end
    ImageComments{counter,1} = currentME.fBIRN_Phantom_QA.SeriesInfo.ImageComments;
    counter = counter+1;
end


%========================================================%

outfile = [parms.outdir '/ABCD_PQA.csv'];
T = table(Site,Manufacturer,ScannerSerialNumber,OperatingSystemVersion,Coil,PhantomCode,StudyDate,StudyTime,SeriesTime,SeriesNumber,SeriesDescription,RepetitionTime,EchoTime,FlipAngle,...,
    Mean,StandardDev,SNR,SFNR,RMS,TemporalDrift,TemporalDrift_per_minute,Max_Absolute_Temporal_Drift,RDC,FWHM_x,FWHM_y,FWHM_z,MeanGhost,...,
    TopGhost,SpatialDrift_PE,ImageFrequency,ImageComments,OSLevel,StudyInstaceUID);

T = sortrows(T,'StudyDate');

writetable(T,outfile,'Delimiter',',','QuoteStrings',true);


return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'indir','/space/syn11/1/data/ABCD/ABCD_QA_proc',[],...
    'outdir','/space/syn11/1/data/ABCD/ABCD_QA_proc/Summary',[],...
    'ProjID','DAL_ABCD_PQA',[],...
    'fname_siteinfo',sprintf('%s/ProjInfo/ABCD/pqa_sites.mat',getenv('HOME')),[],...
  });
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [site, phantom] = lookupSP(st_folder,serial,pqa_sites)
  idx = strcmp({pqa_sites.siteid},st_folder(1:4));
  if ~isempty(serial)
    idx2 = strcmp({pqa_sites.scannersn},serial);
    if ~isempty(find(idx2>0)), idx = idx & idx2; end;
  end
  site = pqa_sites(idx).site;
  phantom = pqa_sites(idx).phantomid;
return;
