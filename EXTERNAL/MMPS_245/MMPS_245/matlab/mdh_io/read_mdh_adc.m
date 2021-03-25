function [adc_data, mdh] = read_mdh_adc(fid,VersionNum)
if (VersionNum==15) | (VersionNum==16)
	[adc_data, mdh] = read_mdh_adc_VA15(fid);
elseif (VersionNum==21)
	[adc_data, mdh] = read_mdh_adc_VA21(fid);
end

