function hdr_szs = read_headerinfo(textfile,rdbm_rev)
%read_headerinfo - Read raw data header information
%
%  hdr_szs = read_headerinfo(textfile,rdbm_rev)
%    textfile - Header information file
%    rdbm_rev - raw header (RDBM) revision number
%    hdr_szs - vector with size of each header section (bytes)
%
%  Assumes the following format
%  #REVNUM
%  headersize
%  offset1
%  offset2
%  offset3
% ...
% #end
% As long as the values between the #REVNUM and #end are the correct values
% and no corruption occurs there, the reader does not care about anything
% else outside

% Copyright (c) 2012 by General Electric Company. All rights reserved.

% Revision History:
% Charles Michelich, 2008-05-11, Updated to support UNIX or DOS end-of-line

% % % textfile
% % % rdbm_rev
fid = fopen(textfile,'r');
% % % fid
flag = 0;
cnt = 1;
hdr_szs(1) = -1;
while ~feof(fid)
    str = fgetl(fid);
    if strncmp(str, '#', 1) && ~strncmp(str, '#end', 4)
        % found a line with #REVNUM ... read it
        if str2double(str(2:length(str))) == rdbm_rev
            % found the rev that I was looking for
            str = str(2:length(str));
            flag = 1;
            cnt = 1;
        end
    elseif strncmp(str, '#end', 4)
        if flag == 1
            %done with reading header offsets
            return;
        end
        %otherwise just continue dont do anything
    end

    if flag == 1
        % Matched requested revision, read sizes
        hdr_szs(cnt) = str2double(str);
        cnt = cnt+1;
    end
end
