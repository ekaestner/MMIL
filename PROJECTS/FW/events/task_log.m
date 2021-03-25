dum = mmil_readtext('/home/halgdev/data/iEEG_NYU/NY158_c/Presentation/NY158_081216_FWIO_logfile.txt','\t');
dum = dum(7:end,3);
for iT = 1:numel(dum)
    if isnumeric(dum{iT})
        dum_num(iT) = 0;
    elseif ~isempty(string_find({dum{iT}(1:2)},{'nw'}))
        dum_num(iT) = 3;
    elseif ~isempty(string_find({dum{iT}(1)},{'O'}))
        dum_num(iT) = 4;
    elseif ~isempty(string_find({dum{iT}},{'npnw'}))
        dum_num(iT) = 5;
    elseif ~isempty(string_find({dum{iT}},{'FF'}))
        dum_num(iT) = 6;
    else
        dum_num(iT) = 0;
    end
end
dum_num(dum_num==0)=[];
dum_num'
