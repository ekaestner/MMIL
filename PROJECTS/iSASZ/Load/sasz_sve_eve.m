function sasz_sve_eve(fcfg,dat)

if ~exist([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/'],'dir'); mkdir([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/']); end

%% Load real events
szv_idn_col = 2;
szv_eve_col = 4;
szv_stm = mmil_readtext([fcfg.clr_fld '/' 'task' '/' 'SZi-actual.csv']);

sza_idn_col = 2;
sza_eve_col = 3;
sza_stm = mmil_readtext([fcfg.clr_fld '/' 'task' '/' 'SAi-actual.csv']);
sza_stm(:,sza_eve_col) = num2cell(cell2mat(sza_stm(:,sza_eve_col))+10);

if strcmpi(fcfg.sbj_nme,'NY226_SA_SZ')
    szv_stm(cell2mat(szv_stm(:,szv_eve_col))==8,:) = [];
end

eve_hld = [sza_stm(:,[sza_idn_col sza_eve_col]) ; szv_stm(:,[szv_idn_col szv_eve_col]) ];

%% Check for missing events
for iL = 1:numel(fcfg.inpath)
    if ~isempty(string_find(fcfg.end_dir(iL),{'eeg'})) || ~isempty(string_find(fcfg.end_dir(iL),{'set'}))
        
        ttt_dat = ts_load_data(fcfg.inpath{iL});
        
        eve_lat{iL} = [];
        eve_epc{iL} = [];
        eve_cde{iL,1} = [];
        
        for iE = 1:numel(ttt_dat.epochs)
            
            eve_lat{iL} = [eve_lat{iL} ttt_dat.epochs(iE).trial_info.latency];
            eve_epc{iL} = [eve_epc{iL} ttt_dat.epochs(iE).trial_info.number];
            eve_cde{iL,1} = [eve_cde{iL,1} ttt_dat.epochs(iE).trial_info.event_code];
            
        end
        
        [~,eve_ind] = sort(eve_epc{iL});
        eve_cde{iL,1} = eve_cde{iL}(eve_ind)';
        eve_lat{iL,1} = eve_lat{iL}(eve_ind)';
       
    elseif  ~isempty(string_find(fcfg.end_dir(iL),{'.mat'}))
        
        nme = strfind(fcfg.inpath{iL},'.');
        nme = fcfg.inpath{iL}(1:nme(1)-1);
        
        fwo = [nme '.dio.txt'];
        try fwo = mmil_readtext(fwo,' ');
        catch
            nme_ind = strfind(nme,'/');
            nme = nme(nme_ind(end)+1:nme_ind(end)+18);
            
            cfg = [];
            cfg.indir   = fcfg.indir;
            cfg.cln_fld = fcfg.cln_fld(iL);
            cfg.cln_dir = {'*'};
            cfg.end_dir = {'dio.txt'};
            [inpath,fwv_dat,trl] = mmil_find_files(cfg);
            fwo = mmil_readtext(inpath{string_find(inpath,{nme})},' ');
        end
        
        if ~any(cell2mat(fwo(1,2:33))); fwo(1,:) = []; end
        fwo = cell2mat(fwo(1:2:end,[1 17:25]));
        
        for iFW = 1:size(fwo,1)
            eve_cde{iL,1}(iFW,1) = bin2dec(num2str(fwo(iFW,2:end)));
            eve_lat{iL,1}(iFW,1) = fwo(iFW,1);
        end
        
        rmv_ind = find(eve_cde{iL,1}==64);
        eve_cde{iL,1}(rmv_ind) = [];
        eve_lat{iL,1}(rmv_ind) = [];        
        
    elseif ~isempty(string_find(fcfg.end_dir(iL),{'edf'}))
        
        try eve_cde{iL,1} = dat.(dat.data_name{iL}).trialinfo; end
        
    end
end

if size(eve_cde,1) ~= numel(fcfg.inpath)
    num_tsk = unique(fcfg.tsk);
    for iN = 1:numel(num_tsk)
        eve_cde{iN,1} = dat.(dat.data_name{iN}).trialinfo;
    end
end

%% Fix event numbers
szv_eve_key = [ 11  10  13  12  35   95   115  54    104  27   47   76   86   66    ; ...
                1   2   3   4   5    5    5    6     6    7    7    8    8    8     ];
sza_eve_key = [ 112 111 122 121 8212 8222 5212 5222 3211 3221  10211 10221 9211 9221 1212 1222 7222 2222 4221 6221 6211 ; ...
                1   2   3   4   5    5    5    5    6    6     6     6     6    6    7    7    7    7    8    8    8 ];

for iL = 1:size(eve_cde,1)
    if any(eve_cde{iL,1}>20)
        if ~isempty(string_find(fcfg.tsk(iL),{'sz'}))
            for iE = 1:size(szv_eve_key,2)
                eve_cde{iL,1}(eve_cde{iL,1}==szv_eve_key(1,iE)) = szv_eve_key(2,iE);
                dat.(dat.data_name{iL}).trialinfo(dat.(dat.data_name{iL}).trialinfo==szv_eve_key(1,iE)) = szv_eve_key(2,iE);
            end
        else
            for iE = 1:size(sza_eve_key,2)
                eve_cde{iL,1}(eve_cde{iL,1}==sza_eve_key(1,iE)) = sza_eve_key(2,iE)+10;
                dat.(dat.data_name{iL}).trialinfo(dat.(dat.data_name{iL}).trialinfo==sza_eve_key(1,iE)) = sza_eve_key(2,iE)+10;
            end
        end
    else
        if ~isempty(string_find(fcfg.tsk(iL),{'sa'}))
            if ~all(eve_cde{iL,1}>10 & eve_cde{iL,1}<20)
                eve_cde{iL,1} = eve_cde{iL,1} + 10;
            dat.(dat.data_name{iL}).trialinfo = dat.(dat.data_name{iL}).trialinfo+10;
            end 
        end
    end
end

%% Save
cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events.csv'],[num2cell(1:size(eve_hld,1))' eve_hld]);

if ~isempty(string_find(fcfg.end_dir(iL),{'eeg'}))  | ~isempty(string_find(fcfg.end_dir(iL),'set'))
    cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes.csv'],[num2cell(cat(1,eve_cde{:})) num2cell(1:numel(num2cell(cat(1,eve_cde{:}))))']);
elseif ~isempty(string_find(fcfg.end_dir(iL),{'mat'}))
    cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes.csv'],[num2cell(cat(1,eve_cde{:})) num2cell(1:numel(num2cell(cat(1,eve_cde{:}))))']);
elseif ~isempty(string_find(fcfg.end_dir(iL),{'edf'}))
    cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_edf_codes.csv'],[num2cell(cat(1,eve_cde{:})) num2cell(1:numel(num2cell(cat(1,eve_cde{:}))))']);
else
    error('')
end

end