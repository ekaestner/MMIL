function [Trans_mat,Cardinal,Extra] = load_trans_apos(fname,verbose)

% function [Trans_mat,Cardinal,Extra] = load_trans_apos(fname,verbose)
%
% LOAD_TRANS_APOS Load in device-to-head transformation information from 
% ".apos" file in fname
% 
% Trans_mat : 4 by 4 device-to-head transformation with rotation and 
% translation 
% Cardinal: coordinates of the LPA NA and RPA in head coord system
% Extra: extra data points in head coord system


fid = fopen(fname,'rt');	% text mode all platforms

if(fid < 0),
  error(['Unable to open file ' fname])
end

next_line = fgetl(fid);
if(verbose)
  disp(next_line)
end
key_word = sscanf(next_line,'%s ',1);

while(~strcmp(key_word,'cardinal')),
  next_line = fgetl(fid);
  if(verbose)
    disp(next_line)
  end
  key_word = sscanf(next_line,'%s ',1);
end				% found keyword CARDINAL

for i = 1:3,
  Cardinal(i,:) = sscanf(next_line,'cardinal %d %f %f %f')';
  next_line = fgetl(fid);
end
if(verbose)
  disp('Read three CARDINAL lines')
end


% so what is the next line?
key_word = sscanf(next_line,'%s ',1);

Extra = [];
while(strcmp(key_word,'extra')), % we have Extra head points
  Extra = [Extra;sscanf(next_line,'extra %d %f %f %f')'];
  next_line = fgetl(fid);
  key_word = sscanf(next_line,'%s ',1);
end
if(verbose),
  disp(sprintf('Read %.0f lines of Extra head points',size(Extra,1)));
end



% now find the transformation matrix

key_word = sscanf(next_line,'%s ',1);
while(strcmp(key_word,'#'))
  next_line = fgetl(fid);
  key_word = sscanf(next_line,'%s ',1);
  if(verbose),
    disp(next_line)
  end
end

for i = 1:3,
  Rot(i,:) = sscanf(next_line,'%f %f %f')';
  next_line = fgetl(fid);
end
Trans = sscanf(next_line,'%f %f %f');

Trans_mat=[[Rot Trans];[0 0 0 1]];

if(verbose),
  disp('Loaded in rotation and translation information')
end


return
