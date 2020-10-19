function outdata = ts_addscript(indata)
%Purpose: Adds date and script fields to timesurfer (or other) data
%         structure. These fields contain character arrays. 
%         Date contains the date ts_addscript was called.
%         Script contains the script (or function) in which ts_addscript was called.
%         Note: ts_addscript may not be used from the command line but must
%         be within an .m file.
%         Note: if ts_addscript is used on a structure already containing
%         the date and/or script fields, the fields are overwritten.
%
% Example: epoch_data = ts_addscript(epoch_data)
%
% Inputs: data structure of any type.
% Outputs: Input data struct with the date and script fields added.
%
%Created by BQR 05/18/12


mfile = dbstack('-completenames');
%% sanity checks
if ~isstruct(indata)
    error('ts_addscript:Input must be struct.')
end
if isempty(mfile)
   error('ts_addscript:ts_addscript may not be used from the command line but must be within an .m file.') 
end
if isfield(indata,'date');
    warning('ts_addscript:date field being overwritten.')
end
if isfield(indata,'script');
    warning('ts_addscript:script field being overwritten.')
end
%% read mfile
mfile = mfile(2).file;
fid = fopen(mfile);
script = fscanf(fid,'%c',Inf);
fclose(fid);
%% print mfile contents to struct field
indata.script = script;
indata.date = date;
outdata = indata;
end