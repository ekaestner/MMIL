clear; clc;



%% Constants
prj_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/iEEG/fce_nme_stm/PreliminaryData/';
out_plt = [ prj_dir '/' 'memory_figures' '/'];

% List names
enc_lst_fld = { 'encoding_associative' };
rcl_lst_fld = { 'recall_associative' };
rec_lst_fld = { 'recall_recognition' };

% Log names
enc_log_fld = { 'encoding_associative' };
rcl_log_fld = { 'recall_associative' };
rec_log_fld = { 'recognition_faces' 'recognition_names' };

% Subjects
sbj_fld = { 'ejk001' };
con_fld = { 'ejk001' };

%% Identify Files
for iS = 1:numel(sbj_fld)
    
    % Encoding Lists
    enc_lst{iS} = dir([ prj_dir '/' sbj_fld{iS} '/' 'lists' '/' enc_lst_fld{1} '/' '*.csv' ]);
    enc_lst{iS} = {enc_lst{iS}(:).name};
    
    % Encoding Log
    enc_log{iS} = dir([ prj_dir '/' sbj_fld{iS} '/' 'logs' '/' enc_log_fld{1} '/' sbj_fld{iS} '*.log' ]);
    enc_log{iS} = {enc_log{iS}(:).name};
    enc_log{iS}(string_find(enc_log{iS},'practice')) = [];
    
    % Recall Lists
    rcl_lst{iS} = dir([ prj_dir '/' sbj_fld{iS} '/' 'lists' '/' rcl_lst_fld{1} '/' '*.csv' ]);
    rcl_lst{iS} = {rcl_lst{iS}(:).name};
    
    % Recall Log
    rcl_log{iS} = dir([ prj_dir '/' sbj_fld{iS} '/' 'logs' '/' rcl_log_fld{1} '/' sbj_fld{iS} '*.log' ]);
    rcl_log{iS} = {rcl_log{iS}(:).name};
    rcl_log{iS}(string_find(rcl_log{iS},'practice')) = [];
    
    % Recognition - Lists
    rec_lst{iS} = dir([ prj_dir '/' sbj_fld{iS} '/' 'lists' '/' rec_lst_fld{1} '/' '*.csv' ]);
    rec_lst{iS} = {rec_lst{iS}(:).name};
    
    % Recognition - Logs
    for iREC = 1:numel(rec_log_fld)
        rec_dta{iS}.(rec_log_fld{iREC}) = dir([ prj_dir '/' sbj_fld{iS} '/' 'logs' '/' rec_log_fld{iREC} '/' sbj_fld{iS} '*.log' ]);
        rec_dta{iS}.(rec_log_fld{iREC}) = {rec_dta{iS}.(rec_log_fld{iREC})(:).name};
        rec_dta{iS}.(rec_log_fld{iREC})(string_find(rec_dta{iS}.(rec_log_fld{iREC}),'practice')) = [];
    end
    
end

%% Calculate d'
for iS = 1:numel(sbj_fld)
    for iREC = 1:numel(rec_log_fld)
        
        if strcmpi(rec_log_fld{iREC},'recognition_faces'); stm_col = 1; elseif strcmpi(rec_log_fld{iREC},'recognition_names'); stm_col = 2; end
        
        dta_hld     = cell(numel(rec_dta{iS}.(rec_log_fld{iREC})),1);
        dta_hld_nme = { 'Stimulus' 'Response' 'TrialType' 'Correct' 'ResponseTime'  };
        
        rec_hld     = cell(numel(rec_dta{iS}.(rec_log_fld{iREC}))+1,7);
        rec_hld_nme = { 'Hit_Rate' 'Hit_RT' 'FalseAlarm_Rate' 'FalseAlarm_RT' 'Rejection_Rate' 'Rejection_RT' 'd''' 'c' };
        
        for iDH = 1:numel(dta_hld)
            % Load List %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            lst_hld = mmil_readtext([ prj_dir '/' sbj_fld{iS} '/' 'lists' '/' rec_lst_fld{1} '/' rec_lst{iS}{iDH}]);
            lst_hld(1,:) = [];
            % Load Log %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            lod_hld = mmil_readtext([ prj_dir '/' sbj_fld{iS} '/' 'logs' '/' rec_log_fld{iREC} '/' rec_dta{iS}.(rec_log_fld{iREC}){iDH} ],['\t']);
            kep_ind = find(strcmpi(lod_hld(:,3),'Picture') | strcmpi(lod_hld(:,3),'Response'));
            lod_hld = lod_hld(kep_ind(2:end),:);
            % Extract %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            pic_trl = find(strcmpi(lod_hld(:,3),'Picture'));
            dta_hld{iDH} = cell(numel(pic_trl),numel(dta_hld_nme));
            for iP = 1:numel(pic_trl)
                % Get stimulus name
                dta_hld{iDH}{iP,1} = lst_hld{iP,stm_col}; % Stimulus
                if strcmpi(lod_hld{pic_trl(iP)+1,3},'Response'); dta_hld{iDH}{iP,2}=lod_hld{pic_trl(iP)+1,4}; else; dta_hld{iDH}{iP,2}=NaN; end % Patient response
                dta_hld{iDH}{iP,3} = lst_hld{iP,17}; % Actual Trial Type
                if strcmpi(lod_hld{pic_trl(iP)+1,3},'Response'); dta_hld{iDH}{iP,5} = (lod_hld{pic_trl(iP)+1,5} - lod_hld{pic_trl(iP),5}) / 10; else; dta_hld{iDH}{iP,5} = NaN; end % Response Time
                % Correct or Not
                if strcmpi(dta_hld{iDH}{iP,3},'OLD') && dta_hld{iDH}{iP,2}==11
                    dta_hld{iDH}{iP,4} = 1;
                elseif strcmpi(dta_hld{iDH}{iP,3},'NOVEL') && dta_hld{iDH}{iP,2}==22
                    dta_hld{iDH}{iP,4} = 1;
                else
                    dta_hld{iDH}{iP,4} = 0;
                end
            end
            % cell2csv(); % Save out
            % Calculate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            prf_hld     = cell2mat(dta_hld{iDH}(:,4));
            rsp_tme_hld = cell2mat(dta_hld{iDH}(:,5));
            rec_hld{iDH,1} = (sum(prf_hld(strcmpi(dta_hld{iDH}(:,3),'OLD'))) / sum(strcmpi(dta_hld{iDH}(:,3),'OLD')))*100; % Hit_Rate
            rec_hld{iDH,2} = nanmean(rsp_tme_hld(strcmpi(dta_hld{iDH}(:,3),'OLD') & prf_hld==1, 1)); % Hit_RT
            rec_hld{iDH,3} = (sum(~(prf_hld(strcmpi(dta_hld{iDH}(:,3),'NOVEL')))) / sum(strcmpi(dta_hld{iDH}(:,3),'NOVEL')))*100;% FalseAlarm_Rate
            rec_hld{iDH,4} = nanmean(rsp_tme_hld(strcmpi(dta_hld{iDH}(:,3),'NOVEL') & prf_hld==0, 1));% FalseAlarm_RT
            rec_hld{iDH,5} = (sum(prf_hld(strcmpi(dta_hld{iDH}(:,3),'NOVEL'))) / sum(strcmpi(dta_hld{iDH}(:,3),'NOVEL')))*100;% Rejection_Rate
            rec_hld{iDH,6} = nanmean(rsp_tme_hld(strcmpi(dta_hld{iDH}(:,3),'NOVEL') & prf_hld==1, 1)); % Rejection_RT
            rec_hld{iDH,7} = norminv((rec_hld{iDH,1}-1)/100)-norminv((rec_hld{iDH,3}+1)/100);% d'
            rec_hld{iDH,8} = -0.5*(norminv((rec_hld{iDH,1}-1)/100)+ norminv((rec_hld{iDH,3}+1)/100));
        end
        
        % PutAllTogether %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Combine extract
        tot_dta_hld = cat(1,dta_hld{:});
        % Combine calculate
        prf_hld     = cell2mat(tot_dta_hld(:,4));
        rsp_tme_hld = cell2mat(tot_dta_hld(:,5));
        rec_hld{end,1} = (sum(prf_hld(strcmpi(tot_dta_hld(:,3),'OLD'))) / sum(strcmpi(tot_dta_hld(:,3),'OLD')))*100; % Hit_Rate
        rec_hld{end,2} = nanmean(rsp_tme_hld(strcmpi(tot_dta_hld(:,3),'OLD') & prf_hld==1, 1)); % Hit_RT
        rec_hld{end,3} = (sum(~(prf_hld(strcmpi(tot_dta_hld(:,3),'NOVEL')))) / sum(strcmpi(tot_dta_hld(:,3),'NOVEL')))*100;% FalseAlarm_Rate
        rec_hld{end,4} = nanmean(rsp_tme_hld(strcmpi(tot_dta_hld(:,3),'NOVEL') & prf_hld==0, 1));% FalseAlarm_RT
        rec_hld{end,5} = (sum(prf_hld(strcmpi(tot_dta_hld(:,3),'NOVEL'))) / sum(strcmpi(tot_dta_hld(:,3),'NOVEL')))*100;% Rejection_Rate
        rec_hld{end,6} = nanmean(rsp_tme_hld(strcmpi(tot_dta_hld(:,3),'NOVEL') & prf_hld==1, 1)); % Rejection_RT
        rec_hld{end,7} = norminv((rec_hld{end,1}-1)/100)-norminv((rec_hld{end,3}+1)/100);% d'
        rec_hld{end,8} = -0.5*(norminv((rec_hld{end,1}-1)/100)+ norminv((rec_hld{end,3}+1)/100));
        % Save out data & list/log match
        cell2csv( [ prj_dir '/' sbj_fld{iS} '/' sbj_fld{iS} '_' rec_log_fld{iREC} '_score.csv' ], rec_hld)
        
    end
end

%% Calculate Recall
for iS = 1:numel(sbj_fld)
    
    if numel(enc_log{iS}) ~= numel(rcl_log{iS}); error('Fix Recall/Encoding Log files'); end
    
    % Fill in File
    for iL = 1:numel(enc_log{iS})
        
        enc_lst_hld{iL} = mmil_readtext( [ prj_dir '/' sbj_fld{iS} '/' 'lists' '/' enc_lst_fld{1} '/' enc_lst{iS}{iL}]);
        kep_ind = find(strcmpi(enc_lst_hld{iL}(:,3),'pl'));
        enc_lst_hld{iL} = enc_lst_hld{iL}(kep_ind,:);
        enc_log_hld{iL} = mmil_readtext( [ prj_dir '/' sbj_fld{iS} '/' 'logs' '/'  enc_log_fld{1} '/' enc_log{iS}{iL} ],['\t'] );
        kep_ind = find(strcmpi(enc_log_hld{iL}(:,4),'255: enc_fce_eve'));
        enc_log_hld{iL} = enc_log_hld{iL}(kep_ind,:);
        rcl_lst_hld{iL} = mmil_readtext( [ prj_dir '/' sbj_fld{iS} '/' 'lists' '/' rcl_lst_fld{1} '/' rcl_lst{iS}{iL}]);
        kep_ind = find(strcmpi(rcl_lst_hld{iL}(:,3),'delayed_recall'));
        rcl_lst_hld{iL} = rcl_lst_hld{iL}(kep_ind,:);
        rcl_log_hld{iL} = mmil_readtext( [ prj_dir '/' sbj_fld{iS} '/' 'logs' '/'  rcl_log_fld{1} '/' rcl_log{iS}{iL} ],['\t'] );
        kep_ind = find(strcmpi(rcl_log_hld{iL}(:,3),'Picture'));
        rcl_log_hld{iL} = rcl_log_hld{iL}(kep_ind(2:end),:);
        
        stm_hld{iL} = cell(size(rcl_log_hld{iL},1),8);
        stm_ind     = 1:size(rcl_log_hld{iL},1):size(enc_log_hld{iL},1);
        
        stm_hld{iL}(:,1) = enc_lst_hld{iL}(stm_ind(1):stm_ind(2)-1,2);
        stm_hld{iL}(:,3) = enc_lst_hld{iL}(stm_ind(2):stm_ind(3)-1,2);
        stm_hld{iL}(:,5) = enc_lst_hld{iL}(stm_ind(3):end,2);
        stm_hld{iL}(:,7) = rcl_lst_hld{iL}(:,2);
        
        stm_hld{iL} = [ stm_hld{iL} ; cell(1,8)];
        
    end
    
    cell2csv([ prj_dir '/' sbj_fld{iS} '/' sbj_fld{iS} '_recall.csv' ],cat(1,stm_hld{:}))
    
end

% Calculate
for iS = 1:numel(sbj_fld)
    
    rcl_beh = cell(1,4);
    
    fll_sht = mmil_readtext([ prj_dir '/' sbj_fld{iS} '/' sbj_fld{iS} '_recall_filled.csv' ]);
    fll_sht(cellfun(@isempty,fll_sht(:,1)),:) = [];
    
    rcl_beh{1,1} = (nansum(cell2mat(fll_sht(:,2))) / size(fll_sht,1)) * 100;
    rcl_beh{1,2} = (nansum(cell2mat(fll_sht(:,4))) / size(fll_sht,1)) * 100;
    rcl_beh{1,3} = (nansum(cell2mat(fll_sht(:,6))) / size(fll_sht,1)) * 100;
    rcl_beh{1,4} = (nansum(cell2mat(fll_sht(:,8))) / size(fll_sht,1)) * 100;
    
    cell2csv([ prj_dir '/' sbj_fld{iS} '/' sbj_fld{iS} '_recall_score.csv' ],rcl_beh)
    
end

%% Plot
% Patients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pat_enc_one = [];
pat_enc_two = [];
pat_enc_thr = [];
pat_rcl_one = [];

pat_hit_fce = [];
pat_fls_fce = [];
pat_hit_nme = [];
pat_fls_nme = [];

pat_hit_rsp_fce = [];
pat_fls_rsp_fce = [];
pat_rej_rsp_fce = [];
pat_hit_rsp_nme = [];
pat_fls_rsp_nme = [];
pat_rej_rsp_nme = [];

pat_dpr_fce = [];
pat_bia_fce = [];
pat_dpr_nme = [];
pat_bia_nme = [];

for iS = 1:numel(sbj_fld)
    rcl_scr     = mmil_readtext([ prj_dir '/' sbj_fld{iS} '/' sbj_fld{iS} '_recall_score.csv' ]);
    rec_nme_scr = mmil_readtext([ prj_dir '/' sbj_fld{iS} '/' sbj_fld{iS} '_' rec_log_fld{2} '_score.csv' ]);
    rec_fce_scr = mmil_readtext([ prj_dir '/' sbj_fld{iS} '/' sbj_fld{iS} '_' rec_log_fld{1} '_score.csv' ]);
    
    pat_enc_one = [ pat_enc_one rcl_scr{1} ];
    pat_enc_two = [ pat_enc_two rcl_scr{2} ];
    pat_enc_thr = [ pat_enc_thr rcl_scr{3} ];
    pat_rcl_one = [ pat_rcl_one rcl_scr{4} ];
    
    pat_hit_fce = [ pat_hit_fce rec_fce_scr{end,1} ];
    pat_fls_fce = [ pat_fls_fce rec_fce_scr{end,3}];
    pat_hit_nme = [ pat_hit_nme rec_nme_scr{end,1} ];
    pat_fls_nme = [ pat_fls_nme rec_nme_scr{end,3} ];
    
    pat_hit_rsp_fce = [ pat_hit_rsp_fce rec_fce_scr{end,2}];
    pat_fls_rsp_fce = [ pat_fls_rsp_fce rec_fce_scr{end,4}];
    pat_rej_rsp_fce = [ pat_rej_rsp_fce rec_fce_scr{end,6}];
    pat_hit_rsp_nme = [ pat_hit_rsp_nme rec_nme_scr{end,2}];
    pat_fls_rsp_nme = [ pat_fls_rsp_nme rec_nme_scr{end,4}];
    pat_rej_rsp_nme = [ pat_rej_rsp_nme rec_nme_scr{end,6}];
    
    pat_dpr_fce = [ pat_dpr_fce rec_fce_scr{end,7}];
    pat_bia_fce = [ pat_bia_fce rec_fce_scr{end,8}];
    pat_dpr_nme = [ pat_dpr_nme rec_nme_scr{end,7}];
    pat_bia_nme = [ pat_bia_nme rec_nme_scr{end,8}];
    
end

% Controls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Recall plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.xdt = { 0.7         1.7         2.7         4.7 };
fcfg.ydt = { pat_enc_one pat_enc_two pat_enc_thr pat_rcl_one };

fcfg.fce_col     = { rgb('light orange') rgb('light orange') rgb('light orange') rgb('light orange')};
fcfg.edg_col     = { [0 0 0] [0 0 0] [0 0 0] [0 0 0]           };
fcfg.box_plt_col = { rgb('dark orange') rgb('dark orange') rgb('dark orange') rgb('dark orange') };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'Encoding 1' 'Encoding 2' 'Encoding 3' '20m Recall'};
fcfg.ylb = {'Recall Accuracy (%)'};
fcfg.ylm = [ 0 100 ];
fcfg.xlm = [ 0 5.5 ];

fcfg.ttl = [ 'Recall Accuracy' ];

fcfg.out_dir = out_plt;
fcfg.out_nme = [ 'Recall_accuracy' ];

ejk_scatter(fcfg)

% Hit/FA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.xdt = { 0.7         1.7         3.7         4.7 };
fcfg.ydt = { pat_hit_fce pat_fls_fce pat_hit_nme pat_fls_nme };

fcfg.fce_col     = { rgb('light orange') rgb('light orange') rgb('light orange') rgb('light orange')};
fcfg.edg_col     = { [0 0 0] [0 0 0] [0 0 0] [0 0 0]           };
fcfg.box_plt_col = { rgb('dark orange') rgb('dark orange') rgb('dark orange') rgb('dark orange') };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'Face Hits' 'Face FA' 'Name Hits' 'Name FA'};
fcfg.ylb = {'20m Recognition Accuracy (%)'};
fcfg.ylm = [ 0 100 ];
fcfg.xlm = [ 0 5.5 ];

fcfg.ttl = [ 'Recognition Accuracy' ];

fcfg.out_dir = out_plt;
fcfg.out_nme = [ 'Recognition_accuracy' ];

ejk_scatter(fcfg)

% d' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.xdt = { 0.7         1.7         3.7         4.7 };
fcfg.ydt = { pat_dpr_fce pat_dpr_nme pat_bia_fce pat_bia_nme };

fcfg.fce_col     = { rgb('light orange') rgb('light orange') rgb('light orange') rgb('light orange')};
fcfg.edg_col     = { [0 0 0] [0 0 0] [0 0 0] [0 0 0]           };
fcfg.box_plt_col = { rgb('dark orange') rgb('dark orange') rgb('dark orange') rgb('dark orange') };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'Face d''' 'Name d''' 'Face C' 'Name C'};
fcfg.ylb = {'20m Recognition Signal Processing'};
fcfg.xlm = [ 0 5.5 ];

fcfg.ttl = [ 'Recognition Signal Processing' ];

fcfg.out_dir = out_plt;
fcfg.out_nme = [ 'Recognition_Signal_Processing' ];

ejk_scatter(fcfg)

% RT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.xdt = { 0.7             1.7             2.7             4.7             5.7             6.7};
fcfg.ydt = { pat_hit_rsp_fce pat_fls_rsp_fce pat_rej_rsp_fce pat_hit_rsp_nme pat_fls_rsp_nme pat_rej_rsp_nme };

fcfg.fce_col     = { rgb('light orange') rgb('light orange') rgb('light orange') rgb('light orange') rgb('light orange') rgb('light orange')};
fcfg.edg_col     = { [0 0 0] [0 0 0] [0 0 0] [0 0 0] [0 0 0] [0 0 0]          };
fcfg.box_plt_col = { rgb('dark orange') rgb('dark orange') rgb('dark orange') rgb('dark orange') rgb('dark orange') rgb('dark orange') };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'Face Hits' 'Face FA' 'Face Rejections' 'Name Hits' 'Name FA' 'Name Rejections'};
fcfg.ylb = {'20m Recognition RT (ms)'};
fcfg.xlm = [ 0 7.5 ];

fcfg.ttl = [ 'Recognition RT' ];

fcfg.out_dir = out_plt;
fcfg.out_nme = [ 'Recognition_RT' ];

ejk_scatter(fcfg)
 
 
 
 
 