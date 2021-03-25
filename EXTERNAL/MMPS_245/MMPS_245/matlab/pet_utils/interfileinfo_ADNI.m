function info = interfileinfo(filename)
%INTERFILEINFO Read metadata from Interfile 3.3 files.
%   INFO = INTERFILEINFO(FILENAME) returns a structure whose fields contain
%   information about images in an Interfile file.  FILENAME is a string
%   that specifies the name of the graphics file.  The file must be in the
%   current directory or in a directory on the MATLAB path.
%   
%   Examples
%   --------
%
%       info = interfileinfo('dyna3.hdr');
%
%   The sample header file in the example above can be found at:
%
%   http://www.keston.com/Phantoms/
%
%   See also INTERFILEREAD.

%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $   $Date: 2005/06/20 03:10:31 $

info = [];
% check header file extension
[fpath,name,ext] = fileparts(filename);
if isempty(ext)
    filename = [filename '.hdr'];
end

% open file for parsing
fid = fopen(filename);
if fid == -1
    err_id = 'Images:interfileinfo:invalidFilename';
    err_msg = sprintf('Can not find header file %s.', filename);
    error(err_id, err_msg);
end

% initialize variables
bad_chars = '!()[]/-_:;';
dates = ['DateOfKeys' 'ProgramDate' 'PatientDob' 'StudyDate','StudyDateDdMmYryr']; % AMD added to deal with screwy ADNI format
times = ['StudyTime' 'ImageStartTime','StudyTimeHhMmSs']; % AMD added to deal with screwy ADNI format
found_header = 0;
found_end = 0;
line_num = 0;

% parse through the file
while (true)
    line_txt = fgetl(fid);
    % stop if no more lines
    if (line_txt == -1) | ~isempty(strfind(line_txt,'!END OF INTERFILE'))
        break;
        
    % skip empty lines
    elseif (sum(isspace(line_txt)) == length(line_txt))
        continue;

    else
        line_num = line_num+1;
        % find index of separator and issue warning if not found
        sep_ind = strfind(line_txt, ':=');
        if (isempty(sep_ind))
            if strcmp(line_txt,'!INTERFILE')
              found_header = 1;
              continue; % Ignore missing := to deal with screwed up ADNI fmt (AMD added 8/20/06) 
            end
            fclose(fid);
            % if no separator on first non-empty line, then not in INTERFILE format
            if isempty(info)
                err_id = 'Images:interfileinfo:invalidFile';
                err_msg = sprintf('%s is not a valid INTERFILE file.', filename);
                error(err_id, err_msg);
                
            % if not on first non-empty line, then invalid expression
            elseif isempty(strfind(line_txt,'!END OF INTERFILE'))
                err_id = 'Images:interfileinfo:noSeparator';
                err_msg = sprintf('Invalid expression in line %s of %s.', num2str(line_num), filename);
                error(err_id, err_msg);
            end
        
        else
            field_str_ind = 1;
            value_str_ind = sep_ind+2;
            field = '';
            
            % parse string to extract field
            while (true)
                [str, count, errmsg, nextindex] = sscanf(line_txt(field_str_ind:sep_ind-1), '%s', 1);
                % check for duplicate header
                if (strcmp(str, '!INTERFILE'))
                    if (found_header == 1)
                        fclose(fid);
                        err_id = 'Images:interfileinfo:duplicateHeader';
                        err_msg = sprintf('Duplicate Interfile header in line %s of %s.', num2str(line_num), filename);
                        error(err_id, err_msg);
                        
                    else
                        found_header = 1;
                    end
                end
                
                % break if no field in rest of string
                if (count == 0)
                    break;
                end
                
                % concatenate strings to form field
                if (strcmp(str, 'ID'))
                    field = [field str];
                    
                else
                    str = lower(str);
                    i = 1;
                    
                    % remove illegal characters
                    while (i <= length(str))
                        k = strfind(bad_chars, str(i));
                        if (~isempty(k))
                            if (k >= 6)
                                str = [str(1:i-1) upper(str(i+1)) str(i+2:length(str))];

                            else
                                str = [str(1:i-1) str(i+1:length(str))];
                            end
                        end

                        i = i+1;
                    end
                    
                    field = [field upper(str(1)) str(2:length(str))];
                end
                
                field_str_ind = field_str_ind+nextindex-1;
            end
            
            % remove extra spaces from beginning of value string
            for i = value_str_ind:length(line_txt)
                if (~isspace(line_txt(i)))
                    break;
                end
            end
            
            value = strcat(line_txt(i:length(line_txt)), '');
            if (strcmp(field, 'VersionOfKeys'))
                if (~strcmp(value, '3.3'))
                    fclose(fid);
                    err_id = 'Images:interfileinfo:unsupportedVersion';
                    err_msg = 'Unsupported version of keys detected.';
                    error(err_id, err_msg);
                end
            end
            
            if isempty(value)
                value = '';
            end
                
            [x, ok] = str2num(value);
            if ((ok ~= 0) && (isempty(strfind(dates, field))) && (isempty(strfind(times, field))))
                value = x;
            end
            
            % close file if end-of-file marker encountered
            if (strcmp(field, 'EndOfInterfile'))
                found_end = 1;
                break;
                
            else
                % check for header
                if (found_header == 0)
                    fclose(fid);
                    err_id = 'Images:interfileinfo:noHeader';
                    err_msg = 'Interfile header not found.';
                    error(err_id, err_msg);

                % store field and value
                elseif (~strcmp(field, 'Interfile'))
                    if (isfield(info, field))
                        if (ischar(info.(field)))
                            info.(field) = {info.(field) value};
                            
                        elseif (iscell(info.(field)))
                            info.(field){length(info.(field))+1} = value;
                            
                        else
                            info.(field) = [info.(field) value];
                        end
                        
                    else
                        info.(field) = value;
                    end
                end
            end
        end
    end
end

% check for end of file marker
%if (found_end == 0)
if 0 % Ignore missing end of file marker (to work with ADNI data) (Added AMD 8/20/06)
    fclose(fid);
    err_id = 'Images:interfileinfo:unexpectedEOF';
    err_msg = 'Unexpected end of file.';
    error(err_id, err_msg);
end

% close file
fclose(fid);
