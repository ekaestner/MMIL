function tme_out = mmil_calc_time(tme_in)
% format = [hour:minute:sec.msec]

tme_in  = strsplit(tme_in,':');
tme_in(3:4)  = strsplit(tme_in{3},'.');
tme_in  = cellfun(@str2num,tme_in);
tme_out = tme_in(1) * 60 * 60 + tme_in(2) * 60 + tme_in(3) + tme_in(4)/1000;

end
