function new_str_hld = num_2_str(num_2_con)

num_2_con = sort(num_2_con);

new_str_hld = '';
clear ttt_beg_hld
ttt_beg_hld(1,1) = num_2_con(1);
ind = 1;
if (numel(num_2_con)-1) < 1
    ttt_beg_hld(ind,1:2) = [num_2_con(1) num_2_con(1)];
elseif (numel(num_2_con)-1) < 2
    iW = 2;
    if num_2_con(iW) - num_2_con(iW-1) ~= 1;
        ind = ind+1;
        ttt_beg_hld(ind,2) = num_2_con(iW);
    else
        ttt_beg_hld(ind,2) = num_2_con(iW);
    end;
else
    for iW = 2:numel(num_2_con)-1;
        if num_2_con(iW) - num_2_con(iW-1) ~= 1
            ttt_beg_hld(ind,2) = num_2_con(iW-1);
            ind = ind+1;
            ttt_beg_hld(ind,1) = num_2_con(iW); 
        end      
    end
    if ttt_beg_hld(end,1) ~= num_2_con(end) && num_2_con(end) - num_2_con(iW) ~= 1
        ttt_beg_hld(end,2) = num_2_con(end-1);
        ttt_beg_hld(end+1,1) = num_2_con(end);
        ttt_beg_hld(end,2) = num_2_con(end);
    else
        ttt_beg_hld(end,2) = num_2_con(end);
    end
end
for iCH = 1:size(ttt_beg_hld,1)
    if any(ttt_beg_hld(iCH,:)==0)
        if ttt_beg_hld(iCH,1) == 0
            ttt_beg_hld(iCH,1) = ttt_beg_hld(iCH,2);
        else
            ttt_beg_hld(iCH,2) = ttt_beg_hld(iCH,1);
        end
    end
end

new_str_hld = '';
new_str_hld = [new_str_hld '['];
for iS = 1:size(ttt_beg_hld,1)
    if ttt_beg_hld(iS,1) == ttt_beg_hld(iS,2)
        new_str_hld = [new_str_hld num2str(ttt_beg_hld(iS,1)) ' '];
    else
        new_str_hld = [new_str_hld num2str(ttt_beg_hld(iS,1)) ':' num2str(ttt_beg_hld(iS,2)) ' '];
    end
end
new_str_hld(end) = ']';

end
