function json = ejk_read_json(file)

% Open the JSON file for reading
fid = fopen(file, 'r');

% Read the entire file into a character array
rawData = fread(fid, inf);
strData = char(rawData');

% Close the file
fclose(fid);

% Decode the JSON string into a MATLAB structure
json = jsondecode(strData);

end