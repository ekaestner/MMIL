function ts_write_stc(filename,vertices,srate,starttime,sol)
% function ts_write_stc(filename,vertices,srate,starttime,sol)
% 
% writes a stc file with the following input parameters:
%
%     filename
%     vertices      vertex indices in a row vector -> see ts_read_dec_file
%     srate         The sampling rate in Hz
%     starttime     The first time point of the epoch in ms (.i.e -500)
%     sol           The solution Matrix [dipoles * time]
%
[fid,message] = fopen(filename,'w','ieee-be');
if fid == -1
    disp(message)
end

% write starttime in ms
fwrite(fid,starttime,'float32');
% write sampling rate in ms
sample_period=1000/srate;
fwrite(fid,sample_period,'float32');
% write number of vertices/sources
fwrite(fid,length(vertices),'uint32');
% write the source vector
for i=1:length(vertices)
    fwrite(fid,vertices(i),'uint32');
end

% write the number of timepts
fwrite(fid,size(sol,2),'uint32');

% write the solution matrix
for time=1:size(sol,2)
    for dipole=1:size(sol,1)
        fwrite(fid,sol(dipole,time),'float32');
    end
end

% close the file
fclose(fid);
