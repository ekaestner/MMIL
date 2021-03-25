function fw_sve_eve(fcfg)

if ~exist([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/'],'dir'); mkdir([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/']); end

%% Load real events
if strcmpi(fcfg.tsk,'fw_1')
    idn_col = 2;
    eve_col = 3;
    eve_fwv = mmil_readtext([fcfg.clr_fld '/' 'task' '/' 'fw_stimuli_1.csv']);
    eve_fwv(cell2mat(eve_fwv(:,3))==7,:) = [];
else
    idn_col = 4;
    eve_col = 3;
    eve_fwv = mmil_readtext([fcfg.clr_fld '/' 'task' '/' 'fw_stimuli_2.csv']);
    eve_fwv(cell2mat(eve_fwv(:,3))==7,:) = [];
end

eve_cde_all = [];

%% Check for missing events
for iL = 1:numel(fcfg.inpath)
    
    if strcmpi(fcfg.end_dir,'eeg') & (strcmpi(fcfg.tsk,'fw_1') | strcmpi(fcfg.tsk,'fw_2'))
        
        ttt_dat = ts_load_data(fcfg.inpath{iL});
        
        eve_lat = [];
        eve_epc = [];
        eve_cde = [];
        
        for iE = 1:numel(ttt_dat.epochs)
            
            eve_lat = [eve_lat ttt_dat.epochs(iE).trial_info.latency];
            eve_epc = [eve_epc ttt_dat.epochs(iE).trial_info.number];
            eve_cde = [eve_cde ttt_dat.epochs(iE).trial_info.event_code];
            
        end
        
        [~,eve_ind] = sort(eve_epc);
        eve_cde = eve_cde(eve_ind)';
        eve_lat = eve_lat(eve_ind)';
        
    elseif strcmpi(fcfg.end_dir,'mat') | strcmpi(fcfg.end_dir,'.mat') & (strcmpi(fcfg.tsk,'fw_2') | strcmpi(fcfg.tsk,'fw_3'))
        
        nme = strfind(fcfg.inpath{iL},'.');
        nme = fcfg.inpath{iL}(1:nme(1)-1);
        
        fwo = [nme '.dio.txt'];
        try fwo = mmil_readtext(fwo,' ');
        catch
            cfg = [];
            cfg.indir   = fcfg.indir;
            cfg.cln_fld = fcfg.cln_fld;
            cfg.cln_dir = {'*'};
            cfg.end_dir = {'dio.txt'};
            try
                [inpath,fwv_dat,trl] = mmil_find_files(cfg);
                fwo = mmil_readtext(inpath{iL},' ');
            catch
                trialf_add      = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'trialf_add');
                fwo = mmil_readtext(trialf_add{iL},' ');
            end
        end
                
        if ~any(cell2mat(fwo(1,2:33))); fwo(1,:) = []; end
        fwo = cell2mat(fwo(1:2:end,[1 17:25]));
        
        for iFW = 1:size(fwo,1)
            eve_cde(iFW,1) = bin2dec(num2str(fwo(iFW,2:end)));
            eve_lat(iFW,1) = fwo(iFW,1);
        end
        
    elseif strcmpi(fcfg.end_dir,'edf') 
        
        load([fcfg.clr_fld '/' 'trialfun_output' '/' fcfg.sbj_nme '.mat']);
        eve_cde = trl_hld{1}(:,4);
        
    end
    
    if ~exist([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme],'dir'); mkdir([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme]); end
    
    if any(eve_cde==4) && any(eve_cde==8) && any(eve_cde==16) && any(eve_cde==32)
        for iFW = 1:size(eve_cde,1)
            eve_cde(iFW) = numel(factor(eve_cde(iFW)))+1;
        end
    end
    
    if sum(eve_cde==3) > sum(eve_cde==8)
        rmv_ind = eve_cde==1 | eve_cde==2 |              eve_cde==7 | eve_cde==48 | eve_cde==64 | eve_cde>128 | eve_cde==68 | eve_cde==69 | eve_cde==70 | eve_cde==72; % eve_cde==0 | 
    elseif sum(eve_cde==3) < sum(eve_cde==8)
        rmv_ind =  eve_cde==1 | eve_cde==2 | eve_cde==3 | eve_cde==7 | eve_cde==48 | eve_cde==64 | eve_cde>128 | eve_cde==68 | eve_cde==69 | eve_cde==70 | eve_cde==72; % eve_cde==0 |
    end
    eve_cde(rmv_ind) = [];
    if sum(eve_cde==3) < sum(eve_cde==8); eve_cde(eve_cde==8) = 3; end
    
    eve_cde_all = [eve_cde_all ; eve_cde];
    
end

if strcmpi(fcfg.end_dir,'eeg') | strcmpi(fcfg.end_dir,'.eeg') | strcmpi(fcfg.end_dir,'edf')
    cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes.csv'],[num2cell(eve_cde_all) num2cell(1:size(eve_cde_all,1))'])
elseif strcmpi(fcfg.end_dir,'mat')  | strcmpi(fcfg.end_dir,'.mat')
    cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes.csv'],[num2cell(eve_cde_all) num2cell(1:size(eve_cde_all,1))'])
end
cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events.csv'],[eve_fwv(:,1:2) num2cell(1:size(eve_fwv,1))' eve_fwv(:,3)])

end