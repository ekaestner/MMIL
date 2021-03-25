function [starttime, srate, vertices, sol]=ts_read_stc(filename)
% function [starttime, srate, vertices, sol]=ts_read_stc(filename)
% 
% reads a stc file with the following output parameters:
%
%     vertices      vertex indices in a row vector -> see ts_read_dec_file
%     srate         The sampling rate in Hz
%     starttime     The first time point of the epoch in ms (.i.e -500)
%     sol           The solution Matrix [dipoles * time]
%
%  Created:  ?
%  Last Mod:  10/02/10 by Don Hagler
%


[fid,message] = fopen(filename,'r','ieee-be');
if fid == -1
  disp(message)
end

STATUS = fseek(fid, 0, 'eof');
filelen= ftell(fid);
STATUS = fseek(fid, 0, 'bof');

% read starttime in ms
[starttime count]=fread(fid,1,'float32');
% read sampling rate in ms
[sample_period count]=fread(fid,1,'float32');
srate = 1000/(sample_period+eps);
% read number of vertices/sources
[vertices_n count]=fread(fid,1,'uint32');
% read the source vector
[vertices count]=fread(fid,vertices_n,'uint32');
% read the number of timepts
[sol_n count]=fread(fid,1,'uint32');

nperdip=(filelen/4-4-vertices_n)/sol_n/vertices_n;

if round(nperdip) ~= nperdip
  fclose(fid);
  error('stc file %s has wrong length',filename);
end

% read the solution matrix
[sol count]=fread(fid,vertices_n*sol_n*nperdip,'float32');
sol=squeeze(reshape(sol,vertices_n,nperdip,sol_n));
% close the file
fclose(fid);
