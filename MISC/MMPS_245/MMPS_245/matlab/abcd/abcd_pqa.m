function abcd_pqa(st_indir, se_outdir)
%function abcd_pqa(st_indir, se_outdir)
%
% Required inputs:
%  st_indir: input directory
%  se_outdir: output directory
%
%  Created:                by Jose Teruel
%  Last Mod by: 09/12/2017 by Feng Xue

se_list = dir(st_indir);
se_list = se_list([se_list.isdir] & ~strncmpi('.', {se_list.name}, 1)); 

fname_error = fullfile(fileparts(se_outdir),'/errorList.txt');
fname = fullfile(fileparts(se_outdir),'/processedList.txt');
fname_currentjson = fullfile(se_outdir,'/QA_metrics.json');

if exist(fname_currentjson, 'file')
    if exist(fname,'file')
        finishedfileID = fopen(fname,'at');
    else
        finishedfileID = fopen(fname,'wt');
    end
    fprintf(finishedfileID,[st_indir,'\n']);
    fclose(finishedfileID);
    return; 
end
if ~isempty(se_list)
    for i=1:length(se_list)
        [error, msg] = processOneSeries(st_indir,se_outdir,se_list(i).name);
        if exist(fname,'file')
            finishedfileID = fopen(fname,'at');
        else
            finishedfileID = fopen(fname,'wt');
        end
        fprintf(finishedfileID,[st_indir,'\n']);
        fclose(finishedfileID);
        
        if error
            if exist(fname_error,'file')
                errorfileID = fopen(fname_error,'at');
            else
                errorfileID = fopen(fname_error,'wt');
            end
            fprintf(errorfileID,[st_indir,'  ->  ', msg, '\n',]);
            fclose(errorfileID);
        end 
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [error, msg] = processOneSeries(st_indir,se_outdir,se_name)
msg ='';
error = 0;
se_indir=fullfile(st_indir,se_name);
dcmfiles = dir(se_indir);
dcmfiles = dcmfiles(~[dcmfiles.isdir] & ~strncmpi('.', {dcmfiles.name}, 1)); 

try
    info = dicominfo(fullfile(se_indir,dcmfiles(1).name));
catch ME
    return;
end

pID = info.PatientID;
if (isfield(info.PatientName, 'FamilyName'));
    pName = info.PatientName.FamilyName;
else
    pName = '';
end
if (isfield(info, 'SeriesDescription'));
    sDes = info.SeriesDescription;
else
    error=1;
    msg = 'SeriesDescription dicom tag was not found';
    return;
end
% Check this is not a human scan or a field map
if (~isempty(strfind(pID, 'HUMAN')) || ~isempty(strfind(pName, 'HUMAN')) || ~isempty(strfind(sDes, 'AP')) ||...,
        ~isempty(strfind(sDes, 'PA')) || ~isempty(strfind(sDes, 'PMU')))
    return;
end
if (~isempty(strfind(info.SeriesDescription, 'fMRI')) || ~isempty(strfind(info.SeriesDescription, 'fBIRN'))...
        || ~isempty(strfind(info.SeriesDescription, 'FBIRN')) || ~isempty(strfind(info.SeriesDescription, 'FMRI')))
    
    roundedTR = round(info.RepetitionTime);
    %Check output directory and creates it if does not exist
    if ~exist(se_outdir, 'dir')
        fprintf('Creating output directory\n');
        [success, message] = mkdir(se_outdir);
        if ~success;
            fprintf('%s -- %s.m:    ERROR: Problem making output directory - %s\n', datestr(now), mfilename, message);
            error = 1;
            msg = 'ERROR: Problem making output directory';
            return;
        end
    end
       
    if roundedTR == 800     
        if (any(strfind(info.Manufacturer, 'GE')) && (length(dcmfiles)<20000))
            fname = fullfile(se_outdir,'/errorlog.txt');
            if exist(fname,'file')
                finishedfileID = fopen(fname,'at');
            else
                finishedfileID = fopen(fname,'wt');
            end
            msg = sprintf('Not enough files for a GE multiband. Number of files: %d \n',length(dcmfiles));
            fprintf(finishedfileID, msg);
            fclose(finishedfileID);
            error = 1;
            return;
        elseif (any(strfind(info.Manufacturer, 'Philips')) && (length(dcmfiles)<20000))
            fname = fullfile(se_outdir,'/errorlog.txt');
            if exist(fname,'file')
                finishedfileID = fopen(fname,'at');
            else
                finishedfileID = fopen(fname,'wt');
            end
            msg = sprintf('Not enough files for a PHILIPS multiband. Number of files: %d \n',length(dcmfiles));
            fprintf(finishedfileID, msg);
            fclose(finishedfileID);
            error = 1;
            return;
        elseif (any(strfind(info.Manufacturer, 'SIEMENS')) && (length(dcmfiles)<500))
            fname = fullfile(se_outdir,'/errorlog.txt');
            if exist(fname,'file')
                finishedfileID = fopen(fname,'at');
            else
                finishedfileID = fopen(fname,'wt');
            end
            msg = sprintf('Not enough files for a SIEMENS multiband. Number of files: %d \n',length(dcmfiles));
            fprintf(finishedfileID, msg);
            fclose(finishedfileID);
            error = 1;
            return;    
        else          
            abcd_pqa_mb_fbirn_qa_app_main(se_indir,se_outdir);      
        end   
        
    elseif roundedTR == 2000
        
        if (any(strfind(info.Manufacturer, 'GE')) && (length(dcmfiles)<6000))
            fname = fullfile(se_outdir,'/errorlog.txt');
            if exist(fname,'file')
                finishedfileID = fopen(fname,'at');
            else
                finishedfileID = fopen(fname,'wt');
            end
            msg = sprintf('Not enough files for a GE non-multiband. Number of files: %d \n',length(dcmfiles));
            fprintf(finishedfileID, msg);
            fclose(finishedfileID);
            error = 1;
            return;
        elseif (any(strfind(info.Manufacturer, 'Philips')) && (length(dcmfiles)<6000))
            fname = fullfile(se_outdir,'/errorlog.txt');
            if exist(fname,'file')
                finishedfileID = fopen(fname,'at');
            else
                finishedfileID = fopen(fname,'wt');
            end
            msg = sprintf('Not enough files for a Philips non-multiband. Number of files: %d \n',length(dcmfiles));
            fprintf(finishedfileID, msg);
            fclose(finishedfileID);
            error = 1;
            return;
        elseif (any(strfind(info.Manufacturer, 'SIEMENS')) && (length(dcmfiles)<200))
            fname = fullfile(se_outdir,'/errorlog.txt');
            if exist(fname,'file')
                finishedfileID = fopen(fname,'at');
            else
                finishedfileID = fopen(fname,'wt');
            end
            msg = sprintf('Not enough files for a SIEMENS non-multiband. Number of files: %d \n',length(dcmfiles));
            fprintf(finishedfileID, msg);
            fclose(finishedfileID);
            error = 1;
            return;    
        else
            abcd_pqa_fbirn_qa_app_main(se_indir,se_outdir);     
   
        end
    
    else
        fname = fullfile(se_outdir,'/errorlog.txt');
        if exist(fname,'file')
            finishedfileID = fopen(fname,'at');
        else
            finishedfileID = fopen(fname,'wt');
        end
        msg = 'fMRI run without TR equal to 800 or 2000. \n';
        fprintf(finishedfileID, msg);
        fclose(finishedfileID);
        error = 1;
        return;
    end
end 

end
