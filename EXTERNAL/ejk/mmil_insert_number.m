function datnew = mmil_insert_number(dat,ind,ins_dat)

datnew = zeros(1,length(dat)+length(ins_dat)) + NaN;
datnew(ind+(0:length(ind)-1)) = ins_dat;
datnew(isnan(datnew)) = dat;

end

