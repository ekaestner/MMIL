function trl_hld = SL_fix_trl(cfg,trl_hld)

if ~exist([cfg.clr_fld '/' 'task' '/' cfg.sbj]); mkdir([cfg.clr_fld '/' 'task' '/' cfg.sbj]); end
if ~exist([cfg.clr_fld '/' 'task' '/' cfg.sbj '/' 'list']); mkdir([cfg.clr_fld '/' 'task' '/' cfg.sbj '/' 'list']); end
if ~exist([cfg.clr_fld '/' 'task' '/' cfg.sbj '/' 'log']);  mkdir([cfg.clr_fld '/' 'task' '/' cfg.sbj '/' 'log']); end

lst_hld = [];
log_hld = [];

blk_hld = {''};
run_hld = {''};

ind = 1;

for iTS = 1:numel(trl_hld)
    
    if strcmpi(cfg.sbj,'NY439_SL')
        blk_ind = regexpi(cfg.inpath{iTS},'\d_\d');
    elseif strcmpi(cfg.sbj,'NY537_SL') || strcmpi(cfg.sbj,'NY540_SL') || strcmpi(cfg.sbj,'NY523_SL') || strcmpi(cfg.sbj,'NY569_SL') || strcmpi(cfg.sbj,'NY590_SL') || strcmpi(cfg.sbj,'NY591_SL') || strcmpi(cfg.sbj,'NY598_SL')
        blk_ind = regexpi(cfg.inpath{iTS},'\d-\d');
    end
    
    blk = num2str(cfg.inpath{iTS}(blk_ind));
    run = num2str(cfg.inpath{iTS}(blk_ind+2));
    
    if ~any(strcmpi(strcat(blk_hld,run_hld),strcat(blk,run)))
        
        blk_hld{ind} = blk;
        run_hld{ind} = run;
        
        % - EJK Check timings & add events here
        switch cfg.sbj
            case 'NY439_SL'
                
                lst_loc = ['/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY439/NY439_SL/SL_6_12_15_NY439/lists' '/' 'SL_Vis1st_' blk '_' run '_1.csv'];
                log_loc = ['/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY439/NY439_SL/SL_6_12_15_NY439/logfiles' '/' 'NY439-SL_Vis1st_' blk '_' run '.log'];
                
                copyfile(lst_loc,[cfg.clr_fld '/' 'task' '/' cfg.sbj '/' 'list' '/' cfg.sbj '_' blk '_' run])
                copyfile(log_loc,[cfg.clr_fld '/' 'task' '/' cfg.sbj '/' 'log' '/' cfg.sbj '_' blk '_' run])
                
                lst = mmil_readtext(lst_loc);
                
                log = mmil_readtext(log_loc,'\t');
                log = log(6:end,:);
                log(strcmpi(log(:,3),'Response'),:) = [];
                log(cellfun(@isstr,log(:,4)),:) = [];
                log(cellfun(@(x) (x==255 | x==0),log(:,4)),:) = [];
                log = cell2mat(log(:,4));
                
            case 'NY537_SL'
                
                lst_loc = ['/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY537/NY537_SL/SL_10_13_15_NY537/lists' '/' 'SL_Vis1st_' blk '_' run '_1.csv'];
                if iTS<=3
                log_loc = ['/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY537/NY537_SL/SL_10_13_15_NY537/logfiles' '/' 'NY537-SL_Vis1st_' blk '_' run '.log'];
                else
                    log_loc = ['/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY537/NY537_SL/SL_10_13_15_NY537/logfiles' '/' 'NY537_2ndday-SL_Vis1st_' blk '_' run '.log'];
                end
                
                copyfile(lst_loc,[cfg.clr_fld '/' 'task' '/' cfg.sbj '/' 'list' '/' cfg.sbj '_' blk '_' run])
                copyfile(log_loc,[cfg.clr_fld '/' 'task' '/' cfg.sbj '/' 'log' '/' cfg.sbj '_' blk '_' run])
                
                lst = mmil_readtext(lst_loc);
                
                log = mmil_readtext(log_loc,'\t');
                log = log(6:end,:);
                log(strcmpi(log(:,3),'Response'),:) = [];
                log(cellfun(@isstr,log(:,4)),:) = [];
                log(cellfun(@(x) (x==255 | x==0),log(:,4)),:) = [];
                log = cell2mat(log(:,4));
                
            case 'NY540_SL'
                
                lst_loc = ['/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY540/NY540_SL/SL_10_13_15_NY540/lists' '/' 'SL_Vis1st_' blk '_' run '_1.csv'];
                log_loc = ['/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY540/NY540_SL/SL_10_13_15_NY540/logfiles' '/' 'NY540-SL_Vis1st_' blk '_' run '.log'];
                
                copyfile(lst_loc,[cfg.clr_fld '/' 'task' '/' cfg.sbj '/' 'list' '/' cfg.sbj '_' blk '_' run])
                copyfile(log_loc,[cfg.clr_fld '/' 'task' '/' cfg.sbj '/' 'log' '/' cfg.sbj '_' blk '_' run])
                
                lst = mmil_readtext(lst_loc);
                
                log = mmil_readtext(log_loc,'\t');
                log = log(6:end,:);
                log(strcmpi(log(:,3),'Response'),:) = [];
                log(cellfun(@isstr,log(:,4)),:) = [];
                log(cellfun(@(x) (x==255 | x==0),log(:,4)),:) = [];
                log = cell2mat(log(:,4));
                
            case 'NY523_SL'
                
                lst = mmil_readtext('/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY523/NY523_SL/SL_2_11_16_NY523/lists/SL_Vis1st_1_1_1.csv');
                lst2 = mmil_readtext('/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY523/NY523_SL/SL_2_11_16_NY523/lists/SL_Vis1st_1_2_1.csv');
                log_loc = ['/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY523/NY523_SL/SL_2_11_16_NY523/logfiles/NY523' '/' 'NY523-SL_Vis1st_1_1.log'];
                log_loc2 = ['/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY523/NY523_SL/SL_2_11_16_NY523/logfiles/NY523' '/' 'NY523_SL_Vis1st_1_2.log'];
                
                copyfile('/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY523/NY523_SL/SL_2_11_16_NY523/lists/SL_Vis1st_1_1_1.csv',[cfg.clr_fld '/' 'task' '/' cfg.sbj '/' 'list' '/' cfg.sbj '_1_1'])
                copyfile('/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY523/NY523_SL/SL_2_11_16_NY523/lists/SL_Vis1st_1_2_1.csv',[cfg.clr_fld '/' 'task' '/' cfg.sbj '/' 'log' '/' cfg.sbj '_1_2'])
                copyfile(log_loc,[cfg.clr_fld '/' 'task' '/' cfg.sbj '/' 'log' '/' cfg.sbj '_1_1'])
                copyfile(log_loc2,[cfg.clr_fld '/' 'task' '/' cfg.sbj '/' 'log' '/' cfg.sbj '_1_2'])
                
                lst = [lst ; lst2];
                
                log = mmil_readtext(log_loc,'\t');
                log = log(6:end,:);
                log(strcmpi(log(:,3),'Response'),:) = [];
                log(cellfun(@isstr,log(:,4)),:) = [];
                log(cellfun(@(x) (x==255 | x==0),log(:,4)),:) = [];
                log = cell2mat(log(:,4));
                
            case 'NY569_SL'
                
                lst_loc = ['/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY569/NY569_SL/SL_2_11_16_NY569_Font60/lists' '/' 'SL_Vis1st_' blk '_' run '_1.csv'];
                log_loc = ['/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY569/NY569_SL/Logfiles' '/' 'NY569-SL_Vis1st_' blk '_' run '.log'];
                
                copyfile(lst_loc,[cfg.clr_fld '/' 'task' '/' cfg.sbj '/' 'list' '/' cfg.sbj '_' blk '_' run])
                copyfile(log_loc,[cfg.clr_fld '/' 'task' '/' cfg.sbj '/' 'log' '/' cfg.sbj '_' blk '_' run])
                
                lst = mmil_readtext(lst_loc);
                
                log = mmil_readtext(log_loc,'\t');
                log = log(6:end,:);
                log(strcmpi(log(:,3),'Response'),:) = [];
                log(cellfun(@isstr,log(:,4)),:) = [];
                log(cellfun(@(x) (x==255 | x==0),log(:,4)),:) = [];
                log = cell2mat(log(:,4));
                
            case 'NY590_SL'
                
                lst_loc = ['/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY590/NY590_SL/SL_2_11_16/lists' '/' 'SL_Vis1st_' blk '_' run '_1.csv'];
                log_loc = ['/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY590/NY590_SL/SL_2_11_16/logfiles' '/' 'NY590-SL_Vis1st_' blk '_' run '.log'];
                
                copyfile(lst_loc,[cfg.clr_fld '/' 'task' '/' cfg.sbj '/' 'list' '/' cfg.sbj '_' blk '_' run])
                copyfile(log_loc,[cfg.clr_fld '/' 'task' '/' cfg.sbj '/' 'log' '/' cfg.sbj '_' blk '_' run])
                
                lst = mmil_readtext(lst_loc);
                
                log = mmil_readtext(log_loc,'\t');
                log = log(6:end,:);
                log(strcmpi(log(:,3),'Response'),:) = [];
                log(cellfun(@isstr,log(:,4)),:) = [];
                log(cellfun(@(x) (x==255 | x==0),log(:,4)),:) = [];
                log = cell2mat(log(:,4));
                
            case 'NY591_SL'
                
                lst_loc = ['/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY591/NY591_SL/SL_2_11_16/lists' '/' 'SL_Vis1st_' blk '_' run '_1.csv'];
                log_loc = ['/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY591/NY591_SL/SL_2_11_16/logfiles/NY591/' '/' 'NY591-SL_Vis1st_' blk '_' run '.log'];
                
                copyfile(lst_loc,[cfg.clr_fld '/' 'task' '/' cfg.sbj '/' 'list' '/' cfg.sbj '_' blk '_' run])
                copyfile(log_loc,[cfg.clr_fld '/' 'task' '/' cfg.sbj '/' 'log' '/' cfg.sbj '_' blk '_' run])
                
                lst = mmil_readtext(lst_loc);
                
                log = mmil_readtext(log_loc,'\t');
                log = log(6:end,:);
                log(strcmpi(log(:,3),'Response'),:) = [];
                log(cellfun(@isstr,log(:,4)),:) = [];
                log(cellfun(@(x) (x==255 | x==0),log(:,4)),:) = [];
                log = cell2mat(log(:,4));
                
            case 'NY598_SL'
                
                lst_loc = ['/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY598/NY598_SL/SL_2_11_16_NY598/lists/' '/' 'SL_Vis1st_' blk '_' run '_1.csv'];
                log_loc = ['/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY598/NY598_SL/SL_2_11_16_NY598/logfiles/' '/' 'NY598-SL_Vis1st_' blk '_' run '.log'];
                
                copyfile(lst_loc,[cfg.clr_fld '/' 'task' '/' cfg.sbj '/' 'list' '/' cfg.sbj '_' blk '_' run])
                copyfile(log_loc,[cfg.clr_fld '/' 'task' '/' cfg.sbj '/' 'log' '/' cfg.sbj '_' blk '_' run])
                
                lst = mmil_readtext(lst_loc);
                
                log = mmil_readtext(log_loc,'\t');
                log = log(6:end,:);
                log(strcmpi(log(:,3),'Response'),:) = [];
                log(cellfun(@isstr,log(:,4)),:) = [];
                log(cellfun(@(x) (x==255 | x==0),log(:,4)),:) = [];
                log = cell2mat(log(:,4));
                
        end
        
        lst_hld = [lst_hld ; lst];
        log_hld = [log_hld ; log];
        
        ind = ind + 1;
        
    end
    
end

cell2csv([cfg.clr_fld '/' 'task' '/' cfg.sbj '/' cfg.sbj '_list.csv'],[num2cell(1:size(lst_hld,1))' lst_hld])
cell2csv([cfg.clr_fld '/' 'task' '/' cfg.sbj '/' cfg.sbj '_log.csv'],num2cell(log_hld))

end