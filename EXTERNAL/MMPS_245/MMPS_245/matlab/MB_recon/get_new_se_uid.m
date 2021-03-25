function [se_uid] = get_new_se_uid(dcmfile);

dicomdict('set','gems-dicom-dict.txt');
metadata = dicominfo(dcmfile);

old_seuid = metadata.SeriesInstanceUID;
splitStr = regexp(old_seuid, '\.', 'split');
addOn = floor(rand(1)*100);
splitStr{end} = num2str(str2double(splitStr{end})+addOn);

se_uid = splitStr{1};
for n = 2:numel(splitStr)
    se_uid = strcat(se_uid,'.',splitStr{n});
end