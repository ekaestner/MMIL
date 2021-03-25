function [Sensor_info,Trans_mat,Cardinal,Extra,EEG_position,HPI_patient,HPI_device] = load_fif_306_apos(fname,verbose)


%function [Sensor_info,Trans_mat,Cardinal,Extra,EEG_position,HPI_patient,HPI_device] = load_fif_306_apos(fname,verbose)
%
% LOAD_FIF_306_APOS Load in position information from Vector-view 306 ".apos" file
% Rcoils is the coil locations, Ocoils is the orientation
% Rot is the rotation matrix, Trans the translation.  If r is a 3 x 1, then
%  r_head = Rot*r + Trans  rotates and translates r from device coordinates
%  to patient coordinates.  The patient coordinates are established by
%  Cardinal, and the HPI coils are given in both the patient and device
%  coordinates.  Sensor_info loads up the device location and three cosine
%  directions for each sensor. Sensor_type loads the type of MEG sensors. 
%  If the last digit in Sensor_type is 1, the sensor is magnetometer, otherwise
%  the sensor is a gradiometer
% 
%  


if(exist('verbose') ~= 1),
  verbose = 1;			% talkative running
end

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

for i = 1:4,
  HPI_patient(i,:) = sscanf(next_line,'hpi %d %f %f %f')';
  next_line = fgetl(fid);
end
if(verbose)
  disp('Read four HPI lines')
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

% so what is the next line?
key_word = sscanf(next_line,'%s ',1);

EEG_position = [];
while(strcmp(key_word,'eeg')), % we have EEG electrodes
  EEG_position = [EEG_position;sscanf(next_line,'eeg %d %f %f %f')'];
  next_line = fgetl(fid);
  key_word = sscanf(next_line,'%s ',1);
end
if(verbose),
  disp(sprintf('Read %.0f lines of EEG electrodes',size(EEG_position,1)));
end


% we are now entering the region where the hpi in device is listed

while(~strcmp(key_word,'hpi'))
  next_line = fgetl(fid);
  key_word = sscanf(next_line,'%s ',1);
  if(verbose),
    disp(next_line)
  end
end

for i = 1:4,
  HPI_device(i,:) = sscanf(next_line,'hpi %d %f %f %f')';
  next_line = fgetl(fid);
end

if(verbose),
  disp('Read in four hpi locations in device coords')
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

% now load in the sensor position information
if(verbose),
  disp('Now loading in sensor information')
end

next_line = fgetl(fid);
key_word = sscanf(next_line,'%s ',inf);

i=0;

while(next_line ~= -1)
      key_word = sscanf(next_line,'%s ',inf);
   if(verbose),
    	disp(next_line)
  	end
   if length(key_word) == 8 | length(key_word) == 7 ,
      if strcmp(key_word(1:4),'#MEG') | strcmp(key_word(1:4),'#EEG'),
         i=i+1;
			Sensor_info(i).name = key_word;
      end   
   elseif strcmp(key_word(1),'#')~=1,
      Sensor_info(i).info = sscanf(next_line,'%f ')';
   end   
   next_line = fgetl(fid);
 end


fclose(fid);

%[Rcoils_grad,Ocoils_grad,Rcoils_mag,Ocoils_mag,x,y,z] = sensor_306(Sensor_info,Sensor_type);

%Rcoils_grad(:,1:3) = [Rot*Rcoils_grad(:,1:3)' + Trans*ones(1,size(Rcoils_grad,1))]';
%Rcoils_grad(:,4:6) = [Rot*Rcoils_grad(:,4:6)' + Trans*ones(1,size(Rcoils_grad,1))]';

%Ocoils_grad(:,1:3) = [Rot*Ocoils_grad(:,1:3)']';
%Ocoils_grad(:,4:6) = [Rot*Ocoils_grad(:,4:6)']';

%Rcoils_mag(:,1:3) = [Rot*Rcoils_mag(:,1:3)' + Trans*ones(1,size(Rcoils_mag,1))]';
%Rcoils_mag(:,4:6) = [Rot*Rcoils_mag(:,4:6)' + Trans*ones(1,size(Rcoils_mag,1))]';
%Rcoils_mag(:,7:9) = [Rot*Rcoils_mag(:,7:9)' + Trans*ones(1,size(Rcoils_mag,1))]';
%Rcoils_mag(:,10:12) = [Rot*Rcoils_mag(:,10:12)' + Trans*ones(1,size(Rcoils_mag,1))]';

%Ocoils_mag(:,1:3) = [Rot*Ocoils_mag(:,1:3)']';
%Ocoils_mag(:,4:6) = [Rot*Ocoils_mag(:,4:6)']';
%Ocoils_mag(:,7:9) = [Rot*Ocoils_mag(:,7:9)']';
%Ocoils_mag(:,10:12) = [Rot*Ocoils_mag(:,10:12)']';

%weight_grad = [1 -1]*16.8/1000;
%weight_mag = [1 1 1 1]/4;

return
