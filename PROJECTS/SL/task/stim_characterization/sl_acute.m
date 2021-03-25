%% Acute
clear; clc;

beh_dir = '/home/ekaestne/Downloads/sd008_behavior';

beh_fle = dir(beh_dir); beh_fle = {beh_fle(~[beh_fle(:).isdir]).name};

for iB = 1:numel(beh_fle)
   
    beh_ana_fle = mmil_readtext([beh_dir '/' beh_fle{iB}],'\t');
    beh_ana_fle = beh_ana_fle(6:end,3:4);
    
    pic_loc = find(strcmpi(beh_ana_fle(:,1),'Picture'));
    
    beh_hld.(mmil_spec_char(beh_fle{iB}(1:end-4),{'-'})).eve_1 = [0 0];
    beh_hld.(mmil_spec_char(beh_fle{iB}(1:end-4),{'-'})).eve_2 = [0 0];
    beh_hld.(mmil_spec_char(beh_fle{iB}(1:end-4),{'-'})).eve_3 = [0 0];
    beh_hld.(mmil_spec_char(beh_fle{iB}(1:end-4),{'-'})).eve_4 = [0 0];
    
    for iT = 1:numel(pic_loc) 
        
        beh_hld.(mmil_spec_char(beh_fle{iB}(1:end-4),{'-'})).(['eve_' num2str(beh_ana_fle{pic_loc(iT),2})])(1) = beh_hld.(mmil_spec_char(beh_fle{iB}(1:end-4),{'-'})).(['eve_' num2str(beh_ana_fle{pic_loc(iT),2})])(1) + 1;
        if beh_ana_fle{pic_loc(iT),2} ~= 1
            if ~strcmpi(beh_ana_fle(pic_loc(iT),1),'Response')
                beh_hld.(mmil_spec_char(beh_fle{iB}(1:end-4),{'-'})).(['eve_' num2str(beh_ana_fle{pic_loc(iT),2})])(2) = beh_hld.(mmil_spec_char(beh_fle{iB}(1:end-4),{'-'})).(['eve_' num2str(beh_ana_fle{pic_loc(iT),2})])(2) + 1;
            end
        else
            if strcmpi(beh_ana_fle(pic_loc(iT)+1,1),'Response')
                beh_hld.(mmil_spec_char(beh_fle{iB}(1:end-4),{'-'})).(['eve_' num2str(beh_ana_fle{pic_loc(iT),2})])(2) = beh_hld.(mmil_spec_char(beh_fle{iB}(1:end-4),{'-'})).(['eve_' num2str(beh_ana_fle{pic_loc(iT),2})])(2) + 1;
            end
        end
    end
    
end

cor = 0;
tot = 0;
beh_nme = fieldnames(beh_hld);
for iC = 1:numel(beh_nme)
    trl_nme = fieldnames(beh_hld.(beh_nme{iC}));
    for iCR = 1:numel(trl_nme)
        cor = cor + beh_hld.(beh_nme{iC}).(trl_nme{iCR})(2);
        tot = tot + beh_hld.(beh_nme{iC}).(trl_nme{iCR})(1);
    end
end
    
for iPL = 1:numel(fieldnames(beh_hld))
    subplot(3,3,iPL)
    bar([beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_1(1)/beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_1(1) ...
         beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_2(1)/beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_2(1) ...
         beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_3(1)/beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_3(1) ...
         beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_4(1)/beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_4(1)],'r');
    hold on;
        bar([beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_1(2)/beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_1(1) ...
             beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_2(2)/beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_2(1) ...
             beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_3(2)/beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_3(1) ...
             beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_4(2)/beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_4(1)],'g');
    title(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-' '_'},'-'))
end

%% Semi-Chronic
clear; clc;

beh_dir = '/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY523/NY523_SL/Logfiles';

beh_fle = dir(beh_dir); beh_fle = {beh_fle(~[beh_fle(:).isdir]).name};
beh_fle = beh_fle(string_find(beh_fle,'_1_'));
beh_fle = beh_fle(string_find(beh_fle,'log'));

for iB = 1:numel(beh_fle)
   
    beh_ana_fle = mmil_readtext([beh_dir '/' beh_fle{iB}],'\t');
    beh_ana_fle = beh_ana_fle(6:end,3:4);
    beh_ana_fle(cellfun(@isstr,beh_ana_fle(:,2)),:) = [];
    
    pic_loc = find(strcmpi(beh_ana_fle(:,1),'Picture'));
    
    beh_hld.(mmil_spec_char(beh_fle{iB}(1:end-4),{'-'})).eve_1 = [0 0];
    beh_hld.(mmil_spec_char(beh_fle{iB}(1:end-4),{'-'})).eve_2 = [0 0];
    beh_hld.(mmil_spec_char(beh_fle{iB}(1:end-4),{'-'})).eve_3 = [0 0];
    beh_hld.(mmil_spec_char(beh_fle{iB}(1:end-4),{'-'})).eve_4 = [0 0];
    
    for iT = 1:numel(pic_loc) 
        
        beh_hld.(mmil_spec_char(beh_fle{iB}(1:end-4),{'-'})).(['eve_' num2str(beh_ana_fle{pic_loc(iT),2})])(1) = beh_hld.(mmil_spec_char(beh_fle{iB}(1:end-4),{'-'})).(['eve_' num2str(beh_ana_fle{pic_loc(iT),2})])(1) + 1;
        if beh_ana_fle{pic_loc(iT),2} ~= 1
            if ~strcmpi(beh_ana_fle(pic_loc(iT),1),'Response')
                if beh_ana_fle{pic_loc(iT)+1,2} == 12
                    beh_hld.(mmil_spec_char(beh_fle{iB}(1:end-4),{'-'})).(['eve_' num2str(beh_ana_fle{pic_loc(iT),2})])(2) = ...
                        beh_hld.(mmil_spec_char(beh_fle{iB}(1:end-4),{'-'})).(['eve_' num2str(beh_ana_fle{pic_loc(iT),2})])(2) + 1;
                    
                end
            end
        else
            if strcmpi(beh_ana_fle(pic_loc(iT)+1,1),'Response')
                if beh_ana_fle{pic_loc(iT)+1,2} == 11
                    beh_hld.(mmil_spec_char(beh_fle{iB}(1:end-4),{'-'})).(['eve_' num2str(beh_ana_fle{pic_loc(iT),2})])(2) = ...
                        beh_hld.(mmil_spec_char(beh_fle{iB}(1:end-4),{'-'})).(['eve_' num2str(beh_ana_fle{pic_loc(iT),2})])(2) + 1;                    
                end
            end
        end
    end
    
end

for iPL = 1:numel(fieldnames(beh_hld))
    subplot(3,3,iPL)
    bar([beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_1(1)/beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_1(1) ...
         beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_2(1)/beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_2(1) ...
         beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_3(1)/beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_3(1) ...
         beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_4(1)/beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_4(1)],'r');
    hold on;
        bar([beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_1(2)/beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_1(1) ...
             beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_2(2)/beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_2(1) ...
             beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_3(2)/beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_3(1) ...
             beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_4(2)/beh_hld.(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-'})).eve_4(1)],'g');
    title(mmil_spec_char(beh_fle{iPL}(1:end-4),{'-' '_'},'-'))
end





















