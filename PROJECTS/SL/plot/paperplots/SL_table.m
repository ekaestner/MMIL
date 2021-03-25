function SL_table

%% Included Regions
inc_reg = {  'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Caudal ITG' 'Middle ITG' 'Rostral ITG' 'Entorhinal' 'Parahippocampal' '' ...
             'Inferior Parietal' 'Superior Parietal' ...
             'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
             'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' 'Supramarginal' ...
             'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

row_cmb = { { 'Lateral Occipital'} ; ...
            { 'Caudal Fusiform' 'Middle Fusiform'} ; ...
            { 'Caudal ITG' 'Middle ITG'} ;...
            { 'Parahippocampal'} ;...
            { 'Entorhinal'} ;...
            {} ;...
            {'Inferior Parietal'} ;...
            {'Superior Parietal'} ;...
            {} ;...
            { 'Inferior Precentral' 'Middle Precentral'} ;...
            { 'Inferior Postcentral' 'Middle Postcentral'} ;...
            {} ;...
            {'Caudal MTG' 'Middle MTG'} ;...
            {'Caudal STG' 'Middle STG'} ;...
            {'Supramarginal'} ;...
            {} ;...
            {'Pars Opercularis'} ;...
            {'Pars Triangularis'} ;...
            {'Pars Orbitalis'} ;...
            {'Caudal Middle Frontal' 'Middle Middle Frontal'} }; 

row_nme = {'Lateral Occipital' ; ...
           'Fusiform' ; ...
           'ITG' ; ...
           'Parahippocampal' ; ...
           'Entorhinal' ; ...
           '' ; ...
           'Inferior Parietal' ; ...
           'Superior Parietal' ; ...
           '' ; ...
           'Precentral' ; ...
           'Postcentral' ; ...
           '' ; ...
           'MTG' ; ...
           'STG' ; ...
           'Supramarginal' ; ...
           '' ; ...
           'Pars Opercularis' ; ...
           'Pars Triangularis' ; ...
           'Pars Orbitalis' ; ...
           'Middle Frontal' };
        
%% Get Data for Selective Electrodes
loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical';
nme = 'pap_anv_1500';

sig_lef{1} = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_lhs_table_plot']); sig_lbl = sig_lef{1}(1,:);
sig_lef{1} = sig_lef{1}(2:end-1,1:6);
sig_lef{1}(:,1) = cellfun(@(x) x(5:end),sig_lef{1}(:,1),'uni',0);

sig_lef{2} = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_lhs_table_plot']); sig_lbl = sig_lef{2}(1,:);
sig_lef{2} = sig_lef{2}(2:end-1,[1 7:11]);
sig_lef{2}(:,1) = cellfun(@(x) x(5:end),sig_lef{2}(:,1),'uni',0);

sig_rgh{1} = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_rhs_table_plot']);
sig_rgh{1} = sig_rgh{1}(2:end-1,1:6);
sig_rgh{1}(:,1) = cellfun(@(x) x(5:end),sig_rgh{1}(:,1),'uni',0);

sig_rgh{2} = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_rhs_table_plot']);
sig_rgh{2} = sig_rgh{2}(2:end-1,[1 7:11]);
sig_rgh{2}(:,1) = cellfun(@(x) x(5:end),sig_rgh{2}(:,1),'uni',0);

% rmv_ind = [];
% for iR = 1:size(sig_lef{1},1)
%     if sig_lef{1}{iR,5}<2
%         rmv_ind = [rmv_ind iR];
%     end
% end
% sig_lef{1}(rmv_ind,:) = [];
% 
% rmv_ind = [];
% for iR = 1:size(sig_lef{2},1)
%     if sig_lef{2}{iR,5}<2
%         rmv_ind = [rmv_ind iR];
%     end
% end
% sig_lef{2}(rmv_ind,:) = [];
% 
% rmv_ind = [];
% for iR = 1:size(sig_rgh{1},1)
%     if sig_rgh{1}{iR,5}<2
%         rmv_ind = [rmv_ind iR];
%     end
% end
% sig_rgh{1}(rmv_ind,:) = [];
% 
% rmv_ind = [];
% for iR = 1:size(sig_rgh{2},1)
%     if sig_rgh{2}{iR,5}<2
%         rmv_ind = [rmv_ind iR];
%     end
% end
% sig_rgh{2}(rmv_ind,:) = [];

% Make Table For Selective Electrodes
tot_reg = unique([sig_lef{1}(:,1) ; sig_lef{2}(:,1) ; sig_rgh{1}(:,1) ;  sig_rgh{2}(:,1)]);
tbl_reg = mmil_order_table(tot_reg);
sig_lef{1} = mmil_order_table(sig_lef{1});
sig_lef{2} = mmil_order_table(sig_lef{2});
sig_rgh{1} = mmil_order_table(sig_rgh{1});
sig_rgh{2} = mmil_order_table(sig_rgh{2});

tbl_reg(strcmpi(tbl_reg,'no_location')) = [];
tbl_reg(strcmpi(tbl_reg,'unknown')) = [];

ovr_tbl_lbl = [repmat({''},1,1) repmat({''},1,1) repmat({''},1,2) {'Left Hemisphere'} repmat({''},1,2) repmat({''},1,1) repmat({''},1,2) {'Right Hemisphere'} repmat({''},1,2)];
tbl_lbl     = {'Region' '' 'Language-Selective Electrodes (Visual 0-450ms)' 'Subjects' '' 'Language-Selective Electrodes (Auditory 450-900ms)' 'Subjects' '' 'Language-Selective Electrodes (Visual 0-450ms)' 'Subjects' '' 'Language-Selective Electrodes (Auditory 450-900ms)' 'Subjects'};
tbl_tbl     = cell(size(tbl_reg,1),13);

for iR = 1:numel(tbl_reg)
   
    tbl_tbl{iR,1} = tbl_reg{iR};
    
    lft_ind = find(strcmpi(sig_lef{1}(:,1),tbl_reg{iR}));
    if ~isempty(lft_ind) > 0
        tbl_tbl{iR,3} = [num2str(round((sig_lef{1}{lft_ind,2}/sig_lef{1}{lft_ind,3})*100)) '% (' [num2str(sig_lef{1}{lft_ind,2}) ' / ' num2str(sig_lef{1}{lft_ind,3})] ')'];
        tbl_tbl{iR,4} = [num2str(sig_lef{1}{lft_ind,5}) ' / ' num2str(sig_lef{1}{lft_ind,6})];
    else
        tbl_tbl{iR,3} = '-';
        tbl_tbl{iR,4} = '-';
    end
    
    lft_ind = find(strcmpi(sig_lef{2}(:,1),tbl_reg{iR}));
    if ~isempty(lft_ind) > 0
        tbl_tbl{iR,6} = [num2str(round((sig_lef{2}{lft_ind,2}/sig_lef{2}{lft_ind,3})*100)) '% (' [num2str(sig_lef{2}{lft_ind,2}) ' / ' num2str(sig_lef{2}{lft_ind,3})] ')'];
        tbl_tbl{iR,7} = [num2str(sig_lef{2}{lft_ind,5}) ' / ' num2str(sig_lef{2}{lft_ind,6})];
    else
        tbl_tbl{iR,6} = '-';
        tbl_tbl{iR,7} = '-';
    end
    
    rgh_ind = find(strcmpi(sig_rgh{1}(:,1),tbl_reg{iR}));
    if ~isempty(rgh_ind) > 0
        tbl_tbl{iR,9} = [num2str(round((sig_rgh{1}{rgh_ind,2}/sig_rgh{1}{rgh_ind,3})*100)) '% (' [num2str(sig_rgh{1}{rgh_ind,2}) ' / ' num2str(sig_rgh{1}{rgh_ind,3})] ')'];
        tbl_tbl{iR,10} = [num2str(sig_rgh{1}{rgh_ind,5}) ' / ' num2str(sig_rgh{1}{rgh_ind,6})];
    else
        tbl_tbl{iR,9} = '-';
        tbl_tbl{iR,10} = '-';
    end
    
    rgh_ind = find(strcmpi(sig_rgh{2}(:,1),tbl_reg{iR}));
    if ~isempty(rgh_ind) > 0
        tbl_tbl{iR,12} = [num2str(round((sig_rgh{2}{rgh_ind,2}/sig_rgh{2}{rgh_ind,3})*100)) '% (' [num2str(sig_rgh{2}{rgh_ind,2}) ' / ' num2str(sig_rgh{2}{rgh_ind,3})] ')'];
        tbl_tbl{iR,13} = [num2str(sig_rgh{2}{rgh_ind,5}) ' / ' num2str(sig_rgh{2}{rgh_ind,6})];
    else
        tbl_tbl{iR,12} = '-';
        tbl_tbl{iR,13} = '-';
    end
    
end

tbl_tbl = mmil_order_table2(tbl_tbl);

rmv_ind = [];
for iTB = 1:size(tbl_tbl,1)
    if ~isempty(tbl_tbl{iTB,1}) && ~ismember(tbl_tbl{iTB,1},inc_reg)
        rmv_ind = [rmv_ind iTB];
    end
end
tbl_tbl(rmv_ind,:) = [];

mmil_chk_dir('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/manuscript/tables')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical';
nme = 'pap_anv_1500';

sig_lef{1} = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_lhs_table_plot']); sig_lbl = sig_lef{1}(1,:);
sig_lef{1} = sig_lef{1}(2:end-1,1:6);
sig_lef{1}(:,1) = cellfun(@(x) x(5:end),sig_lef{1}(:,1),'uni',0);

sig_lef_nme = load([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_lhs_table_plot_subject_names.mat']);
sig_lef_nme.ttt_sav = sig_lef_nme.ttt_sav(2:end-1,1:6);
sig_lef_nme.ttt_sav(1:end,1) = cellfun(@(x) x(5:end),sig_lef_nme.ttt_sav(1:end,1),'uni',0);
sig_lef_nme_hld{1} = sig_lef_nme.ttt_sav;

sig_lef{2} = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_lhs_table_plot']); sig_lbl = sig_lef{2}(1,:);
sig_lef{2} = sig_lef{2}(2:end-1,[1 7:11]);
sig_lef{2}(:,1) = cellfun(@(x) x(5:end),sig_lef{2}(:,1),'uni',0);

sig_lef_nme = load([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_lhs_table_plot_subject_names.mat']);
sig_lef_nme.ttt_sav = sig_lef_nme.ttt_sav(2:end-1,[1 7:11]);
sig_lef_nme.ttt_sav(1:end,1) = cellfun(@(x) x(5:end),sig_lef_nme.ttt_sav(1:end,1),'uni',0);
sig_lef_nme_hld{2} = sig_lef_nme.ttt_sav;

sig_rgh{1} = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_rhs_table_plot']);
sig_rgh{1} = sig_rgh{1}(2:end-1,1:6);
sig_rgh{1}(:,1) = cellfun(@(x) x(5:end),sig_rgh{1}(:,1),'uni',0);

sig_rgh_nme = load([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_rhs_table_plot_subject_names.mat']);
sig_rgh_nme.ttt_sav = sig_rgh_nme.ttt_sav(2:end-1,1:6);
sig_rgh_nme.ttt_sav(1:end,1) = cellfun(@(x) x(5:end),sig_rgh_nme.ttt_sav(1:end,1),'uni',0);
sig_rgh_nme_hld{1} = sig_rgh_nme.ttt_sav;

sig_rgh{2} = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_rhs_table_plot']);
sig_rgh{2} = sig_rgh{2}(2:end-1,[1 7:11]);
sig_rgh{2}(:,1) = cellfun(@(x) x(5:end),sig_rgh{2}(:,1),'uni',0);

sig_rgh_nme = load([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_rhs_table_plot_subject_names.mat']);
sig_rgh_nme.ttt_sav = sig_rgh_nme.ttt_sav(2:end-1,[1 7:11]);
sig_rgh_nme.ttt_sav(1:end,1) = cellfun(@(x) x(5:end),sig_rgh_nme.ttt_sav(1:end,1),'uni',0);
sig_rgh_nme_hld{2} = sig_rgh_nme.ttt_sav;

% rmv_ind = [];
% for iR = 1:size(sig_lef{1},1)
%     if sig_lef{1}{iR,5}<2
%         rmv_ind = [rmv_ind iR];
%     end
% end
% sig_lef{1}(rmv_ind,:) = [];
% 
% rmv_ind = [];
% for iR = 1:size(sig_lef{2},1)
%     if sig_lef{2}{iR,5}<2
%         rmv_ind = [rmv_ind iR];
%     end
% end
% sig_lef{2}(rmv_ind,:) = [];
% 
% rmv_ind = [];
% for iR = 1:size(sig_rgh{1},1)
%     if sig_rgh{1}{iR,5}<2
%         rmv_ind = [rmv_ind iR];
%     end
% end
% sig_rgh{1}(rmv_ind,:) = [];
% 
% rmv_ind = [];
% for iR = 1:size(sig_rgh{2},1)
%     if sig_rgh{2}{iR,5}<2
%         rmv_ind = [rmv_ind iR];
%     end
% end
% sig_rgh{2}(rmv_ind,:) = [];

% Make Table For Selective Electrodes
tot_reg = unique([sig_lef{1}(:,1) ; sig_lef{2}(:,1) ; sig_rgh{1}(:,1) ;  sig_rgh{2}(:,1)]);
tbl_reg = mmil_order_table(tot_reg);
sig_lef{1} = mmil_order_table(sig_lef{1}); sig_lef{1} = mmil_order_table(sig_lef{1});
sig_lef{2} = mmil_order_table(sig_lef{2}); sig_lef{2} = mmil_order_table(sig_lef{2});
sig_rgh{1} = mmil_order_table(sig_rgh{1}); sig_rgh{1} = mmil_order_table(sig_rgh{1});
sig_rgh{2} = mmil_order_table(sig_rgh{2}); sig_rgh{2} = mmil_order_table(sig_rgh{2});

tbl_reg(strcmpi(tbl_reg,'no_location')) = [];
tbl_reg(strcmpi(tbl_reg,'unknown')) = [];

ovr_tbl_lbl = [repmat({''},1,1) repmat({''},1,1) repmat({''},1,2) {'Left Hemisphere'} repmat({''},1,2) repmat({''},1,1) repmat({''},1,2) {'Right Hemisphere'} repmat({''},1,2)];
tbl_lbl     = {'Region' '' 'Language-Selective Electrodes (Visual 0-450ms)' 'Subjects' '' 'Language-Selective Electrodes (Auditory 450-900ms)' 'Subjects' '' 'Language-Selective Electrodes (Visual 0-450ms)' 'Subjects' '' 'Language-Selective Electrodes (Auditory 450-900ms)' 'Subjects'};
tbl_tbl_cmb     = cell(size(row_cmb,1),13);

for iR = 1:numel(row_cmb)
    if ~isempty(row_cmb{iR})
        
    tbl_tbl_cmb{iR,1} = row_nme{iR};
    
    lft_ind     = cellfun(@(x) find(strcmpi(sig_lef{1}(:,1),x)),row_cmb{iR},'uni',0); lft_ind = [lft_ind{:}];
    lft_ind_nme = cellfun(@(x) find(strcmpi(sig_lef_nme_hld{1}(:,1),x)),row_cmb{iR});
    if ~isempty(lft_ind) > 0
        tbl_tbl_cmb{iR,3} = [num2str(round( sum([sig_lef{1}{lft_ind,2}]) / sum([sig_lef{1}{lft_ind,3}]) * 100)) '% (' [num2str( sum([sig_lef{1}{lft_ind,2}]) ) ' / ' num2str( sum([sig_lef{1}{lft_ind,3}]) )] ')'];
        rsp_nme = unique(cat(1,sig_lef_nme_hld{1}{lft_ind_nme,5}));
        tot_nme = unique(cat(1,sig_lef_nme_hld{1}{lft_ind_nme,6}));
        tbl_tbl_cmb{iR,4} = [num2str(numel(rsp_nme)) ' / ' num2str(numel(tot_nme))];
    else
        tbl_tbl_cmb{iR,3} = '-';
        tbl_tbl_cmb{iR,4} = '-';
    end
        
    lft_ind = cellfun(@(x) find(strcmpi(sig_lef{2}(:,1),x)),row_cmb{iR},'uni',0); lft_ind = [lft_ind{:}];
    lft_ind_nme = cellfun(@(x) find(strcmpi(sig_lef_nme_hld{2}(:,1),x)),row_cmb{iR});
    if ~isempty(lft_ind) > 0
        tbl_tbl_cmb{iR,6} = [num2str(round( sum([sig_lef{2}{lft_ind,2}]) / sum([sig_lef{2}{lft_ind,3}]) *100)) '% (' [num2str( sum([sig_lef{2}{lft_ind,2}]) ) ' / ' num2str( sum([sig_lef{2}{lft_ind,3}]) )] ')'];
        rsp_nme = unique(cat(1,sig_lef_nme_hld{2}{lft_ind_nme,5}));
        tot_nme = unique(cat(1,sig_lef_nme_hld{2}{lft_ind_nme,6}));
        tbl_tbl_cmb{iR,7} = [num2str(numel(rsp_nme)) ' / ' num2str(numel(tot_nme))];
    else
        tbl_tbl_cmb{iR,6} = '-';
        tbl_tbl_cmb{iR,7} = '-';
    end
        
    rgh_ind = cellfun(@(x) find(strcmpi(sig_rgh{1}(:,1),x)),row_cmb{iR},'uni',0); rgh_ind = rgh_ind{:};
    if ~isempty(rgh_ind) > 0
        tbl_tbl_cmb{iR,9} = [num2str(round(( sum([sig_rgh{1}{rgh_ind,2}]) / sum([sig_rgh{1}{rgh_ind,3}]) )*100)) '% (' [num2str(sum([sig_rgh{1}{rgh_ind,2}]) ) ' / ' num2str( sum([sig_rgh{1}{rgh_ind,3}]) )] ')'];
        rsp_nme = unique(cat(1,sig_rgh_nme_hld{1}{lft_ind_nme,5}));
        tot_nme = unique(cat(1,sig_rgh_nme_hld{1}{lft_ind_nme,6}));
        tbl_tbl_cmb{iR,10} = [num2str(numel(rsp_nme)) ' / ' num2str(numel(tot_nme))];
    else
        tbl_tbl_cmb{iR,9} = '-';
        tbl_tbl_cmb{iR,10} = '-';
    end
        
    rgh_ind = cellfun(@(x) find(strcmpi(sig_rgh{2}(:,1),x)),row_cmb{iR},'uni',0); rgh_ind = rgh_ind{:};
    if ~isempty(rgh_ind) > 0
        tbl_tbl_cmb{iR,12} = [num2str(round(( sum([sig_rgh{2}{rgh_ind,2}]) / sum([sig_rgh{2}{rgh_ind,3}]) )*100)) '% (' [num2str( sum([sig_rgh{2}{rgh_ind,2}]) ) ' / ' num2str( sum([sig_rgh{2}{rgh_ind,3}]) )] ')'];
        rsp_nme = unique(cat(1,sig_rgh_nme_hld{2}{lft_ind_nme,5}));
        tot_nme = unique(cat(1,sig_rgh_nme_hld{2}{lft_ind_nme,6}));
        tbl_tbl_cmb{iR,13} = [num2str(numel(rsp_nme)) ' / ' num2str(numel(tot_nme))];
    else
        tbl_tbl_cmb{iR,12} = '-';
        tbl_tbl_cmb{iR,13} = '-';
    end
    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmil_chk_dir('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/manuscript/tables')

cell2csv('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/manuscript/tables/table1.csv',[ovr_tbl_lbl ; tbl_lbl ; tbl_tbl])
cell2csv('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/manuscript/tables/table1_cmb.csv',[ovr_tbl_lbl ; tbl_lbl ; tbl_tbl_cmb])

%% Effects Table
loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical';
nme = { 'pap_anv_1500' 'pap_lng_950' 'pap_con_950' 'pap_mtc_1450' 'pap_phn_950'};

% Load Left
anv_lef = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{1} '/' 'total' '/' 'total_' nme{1} '_lhs_table_plot']); sig_lbl = anv_lef(1,:);
anv_lef = anv_lef(2:end-1,[1 17:21]);
anv_lef(:,1) = cellfun(@(x) x(5:end),anv_lef(:,1),'uni',0);

vis_lng_lef = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{2} '/' 'total' '/' 'total_' nme{2} '_lhs_table_plot']); sig_lbl = vis_lng_lef(1,:);
vis_lng_lef = vis_lng_lef(2:end-1,1:6);
vis_lng_lef(:,1) = cellfun(@(x) x(5:end),vis_lng_lef(:,1),'uni',0);

aud_lng_lef = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{2} '/' 'total' '/' 'total_' nme{2} '_lhs_table_plot']); sig_lbl = aud_lng_lef(1,:);
aud_lng_lef = aud_lng_lef(2:end-1,[1 7:11]);
aud_lng_lef(:,1) = cellfun(@(x) x(5:end),aud_lng_lef(:,1),'uni',0);

aud_con_lef = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{3} '/' 'total' '/' 'total_' nme{3} '_lhs_table_plot']); sig_lbl = aud_con_lef(1,:);
aud_con_lef = aud_con_lef(2:end-1,[1 7:11]);
aud_con_lef(:,1) = cellfun(@(x) x(5:end),aud_con_lef(:,1),'uni',0);

vis_con_lef = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{3} '/' 'total' '/' 'total_' nme{3} '_lhs_table_plot']); sig_lbl = vis_con_lef(1,:);
vis_con_lef = vis_con_lef(2:end-1,1:6);
vis_con_lef(:,1) = cellfun(@(x) x(5:end),vis_con_lef(:,1),'uni',0);

inc_lef = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{4} '/' 'total' '/' 'total_' nme{4} '_lhs_table_plot']); sig_lbl = inc_lef(1,:);
inc_lef = inc_lef(2:end-1,[1:6]);
inc_lef(:,1) = cellfun(@(x) x(5:end),inc_lef(:,1),'uni',0);

vis_ltr_lef = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{5} '/' 'total' '/' 'total_' nme{5} '_lhs_table_plot']); sig_lbl = vis_ltr_lef(1,:);
vis_ltr_lef = vis_ltr_lef(2:end-1,1:6);
vis_ltr_lef(:,1) = cellfun(@(x) x(5:end),vis_ltr_lef(:,1),'uni',0);

aud_phn_lef = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{5} '/' 'total' '/' 'total_' nme{5} '_lhs_table_plot']); sig_lbl = aud_phn_lef(1,:);
aud_phn_lef = aud_phn_lef(2:end-1,[1 7:11]);
aud_phn_lef(:,1) = cellfun(@(x) x(5:end),aud_phn_lef(:,1),'uni',0);

% Load Right
anv_rgh = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{1} '/' 'total' '/' 'total_' nme{1} '_rhs_table_plot']); sig_lbl = anv_rgh(1,:);
anv_rgh = anv_rgh(2:end-1,[1 17:21]);
anv_rgh(:,1) = cellfun(@(x) x(5:end),anv_rgh(:,1),'uni',0);

vis_lng_rgh = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{2} '/' 'total' '/' 'total_' nme{2} '_rhs_table_plot']); sig_lbl = vis_lng_rgh(1,:);
vis_lng_rgh = vis_lng_rgh(2:end-1,1:6);
vis_lng_rgh(:,1) = cellfun(@(x) x(5:end),vis_lng_rgh(:,1),'uni',0);

aud_lng_rgh = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{2} '/' 'total' '/' 'total_' nme{2} '_rhs_table_plot']); sig_lbl = aud_lng_rgh(1,:);
aud_lng_rgh = aud_lng_rgh(2:end-1,[1 7:11]);
aud_lng_rgh(:,1) = cellfun(@(x) x(5:end),aud_lng_rgh(:,1),'uni',0);

vis_con_rgh = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{3} '/' 'total' '/' 'total_' nme{3} '_rhs_table_plot']); sig_lbl = vis_con_rgh(1,:);
vis_con_rgh = vis_con_rgh(2:end-1,1:6);
vis_con_rgh(:,1) = cellfun(@(x) x(5:end),vis_con_rgh(:,1),'uni',0);

aud_con_rgh = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{3} '/' 'total' '/' 'total_' nme{3} '_rhs_table_plot']); sig_lbl = aud_con_rgh(1,:);
aud_con_rgh = aud_con_rgh(2:end-1,[1 7:11]);
aud_con_rgh(:,1) = cellfun(@(x) x(5:end),aud_con_rgh(:,1),'uni',0);

inc_rgh = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{4} '/' 'total' '/' 'total_' nme{4} '_rhs_table_plot']); sig_lbl = inc_rgh(1,:);
inc_rgh = inc_rgh(2:end-1,[1:6]);
inc_rgh(:,1) = cellfun(@(x) x(5:end),inc_rgh(:,1),'uni',0);

vis_ltr_rgh = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{5} '/' 'total' '/' 'total_' nme{5} '_rhs_table_plot']); sig_lbl = vis_ltr_rgh(1,:);
vis_ltr_rgh = vis_ltr_rgh(2:end-1,1:6);
vis_ltr_rgh(:,1) = cellfun(@(x) x(5:end),vis_ltr_rgh(:,1),'uni',0);

aud_phn_rgh = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{5} '/' 'total' '/' 'total_' nme{5} '_rhs_table_plot']); sig_lbl = aud_phn_rgh(1,:);
aud_phn_rgh = aud_phn_rgh(2:end-1,[1 7:11]);
aud_phn_rgh(:,1) = cellfun(@(x) x(5:end),aud_phn_rgh(:,1),'uni',0);

% Left Remove
rmv_ind = [];
for iR = 1:size(anv_lef,1)
    if anv_lef{iR,5}<2
        rmv_ind = [rmv_ind iR];
    end
end

% anv_lef(rmv_ind,:) = [];
anv_lef = mmil_order_table(anv_lef);

% vis_lng_lef(rmv_ind,:) = [];
vis_lng_lef = mmil_order_table(vis_lng_lef);

% aud_lng_lef(rmv_ind,:) = [];
aud_lng_lef = mmil_order_table(aud_lng_lef);

% vis_con_lef(rmv_ind,:) = [];
vis_con_lef = mmil_order_table(vis_con_lef);

% aud_con_lef(rmv_ind,:) = [];
aud_con_lef = mmil_order_table(aud_con_lef);

% inc_lef(rmv_ind,:) = [];
inc_lef = mmil_order_table(inc_lef);

vis_ltr_lef = mmil_order_table(vis_ltr_lef); vis_ltr_lef = mmil_order_table(vis_ltr_lef);
aud_phn_lef = mmil_order_table(aud_phn_lef); aud_phn_lef = mmil_order_table(aud_phn_lef);

% Right Remove
rmv_ind = [];
for iR = 1:size(anv_rgh,1)
    if anv_rgh{iR,5}<2
        rmv_ind = [rmv_ind iR];
    end
end

% anv_rgh(rmv_ind,:) = [];
anv_rgh = mmil_order_table(anv_rgh);

% vis_lng_rgh(rmv_ind,:) = [];
vis_lng_rgh = mmil_order_table(vis_lng_rgh);

% aud_lng_rgh(rmv_ind,:) = [];
aud_lng_rgh = mmil_order_table(aud_lng_rgh);

% vis_con_rgh(rmv_ind,:) = [];
vis_con_rgh = mmil_order_table(vis_con_rgh);

% aud_con_rgh(rmv_ind,:) = [];
aud_con_rgh = mmil_order_table(aud_con_rgh);

% inc_rgh(rmv_ind,:) = [];
inc_rgh = mmil_order_table(inc_rgh);

vis_ltr_rgh = mmil_order_table(vis_ltr_rgh); vis_ltr_rgh = mmil_order_table(vis_ltr_rgh);
aud_phn_rgh = mmil_order_table(aud_phn_rgh); aud_phn_rgh = mmil_order_table(aud_phn_rgh);

% Make Table For Selective Electrodes
tot_reg = unique([anv_lef(:,1) ; anv_rgh(:,1)]);
if ~isempty(string_find(tot_reg,{' '})); tbl_reg = tot_reg; else tbl_reg = mmil_order_table(tot_reg); end
tbl_reg(strcmpi(tbl_reg,'no_location')) = [];
tbl_reg(strcmpi(tbl_reg,'unknown')) = [];

ovr_tbl_lbl = [repmat({''},1,2) repmat({''},1,5) {'Left Hemisphere'} repmat({''},1,5) repmat({''},1,1) repmat({''},1,5) {'Right Hemisphere'} repmat({''},1,5)];
tbl_lbl     = {'Region' '' 'Selective Electrodes' '' 'Text Selective' 'False-Font Selective' 'Letter-Sensitive' '' 'Voice Selective' 'Noise-Vocoded Selective' 'Phoneme-Sensitive' '' 'Incongruent Effects' '' 'Selective Electrodes' '' 'Text Selective' 'False-Font Selective' 'Letter-Sensitive' '' 'Voice Selective' 'Noise-Vocoded Selective' 'Phoneme-Sensitive' '' 'Incongruent Effects'};
tbl_tbl     = cell(size(tbl_reg,1),21);

%
for iR = 1:numel(tbl_reg)
    
    tbl_tbl{iR,1} = tbl_reg{iR};
    
    % Left
    lft_ind = find(strcmpi(anv_lef(:,1),tbl_reg{iR}));
    if ~isempty(lft_ind) > 0
        
        lft_ind = find(strcmpi(anv_lef(:,1),tbl_reg{iR}));
        if ~isempty(lft_ind)
            tbl_tbl{iR,3} = anv_lef{lft_ind,2};
        else
            tbl_tbl{iR,3} = 0;
        end
        
        lft_ind = find(strcmpi(vis_lng_lef(:,1),tbl_reg{iR}));
        if ~isempty(lft_ind)
            tbl_tbl{iR,5} = vis_lng_lef{lft_ind,2};
        else
            tbl_tbl{iR,5} = 0;
        end
        
        lft_ind = find(strcmpi(vis_con_lef(:,1),tbl_reg{iR}));
        if ~isempty(lft_ind)
            tbl_tbl{iR,6} = vis_con_lef{lft_ind,2};
        else
            tbl_tbl{iR,6} = 0;
        end
        
        lft_ind = find(strcmpi(vis_ltr_lef(:,1),tbl_reg{iR}));
        if ~isempty(lft_ind)
            tbl_tbl{iR,7} = vis_ltr_lef{lft_ind,2};
        else
            tbl_tbl{iR,7} = 0;
        end
        
        
        lft_ind = find(strcmpi(aud_lng_lef(:,1),tbl_reg{iR}));
        if ~isempty(lft_ind)
            tbl_tbl{iR,9}  = aud_lng_lef{lft_ind,2};
        else
            tbl_tbl{iR,9} = 0;
        end
        
        
        lft_ind = find(strcmpi(aud_con_lef(:,1),tbl_reg{iR}));
        if ~isempty(lft_ind)
            tbl_tbl{iR,10} = aud_con_lef{lft_ind,2};
        else
            tbl_tbl{iR,10} = 0;
        end
        
        
        lft_ind = find(strcmpi(aud_phn_lef(:,1),tbl_reg{iR}));
        if ~isempty(lft_ind)
            tbl_tbl{iR,11} = aud_phn_lef{lft_ind,2};
        else
            tbl_tbl{iR,11} = 0;
        end
        
        lft_ind = find(strcmpi(inc_lef(:,1),tbl_reg{iR}));
        if ~isempty(lft_ind)
            tbl_tbl{iR,13} = inc_lef{lft_ind,2};
        else
            tbl_tbl{iR,13} = 0;
        end
        
    else
        tbl_tbl{iR,3} = '-';
        
        tbl_tbl{iR,5} = '-';
        tbl_tbl{iR,6} = '-';
        tbl_tbl{iR,7} = '-';
        
        tbl_tbl{iR,9}  = '-';
        tbl_tbl{iR,10} = '-';
        tbl_tbl{iR,11} = '-';
        
        tbl_tbl{iR,13} = '-';
    end
    
    % Right
    rgh_ind = find(strcmpi(anv_rgh(:,1),tbl_reg{iR}));
    if ~isempty(rgh_ind) > 0
        
        rgh_ind = find(strcmpi(anv_rgh(:,1),tbl_reg{iR}));
        if ~isempty(rgh_ind)
            tbl_tbl{iR,15} = anv_rgh{rgh_ind,2};
        else
            tbl_tbl{iR,15} = 0;
        end
        
        rgh_ind = find(strcmpi(vis_lng_rgh(:,1),tbl_reg{iR}));
        if ~isempty(rgh_ind)
            tbl_tbl{iR,17} = vis_lng_rgh{rgh_ind,2};
        else
            tbl_tbl{iR,17} = 0;
        end
        
        rgh_ind = find(strcmpi(vis_con_rgh(:,1),tbl_reg{iR}));
        if ~isempty(rgh_ind)
            tbl_tbl{iR,18} = vis_con_rgh{rgh_ind,2};
        else
            tbl_tbl{iR,18} = 0;
        end
        
        rgh_ind = find(strcmpi(vis_ltr_rgh(:,1),tbl_reg{iR}));
        if ~isempty(rgh_ind)
            tbl_tbl{iR,19} = vis_ltr_rgh{rgh_ind,2};
        else
            tbl_tbl{iR,19} = 0;
        end
        
        rgh_ind = find(strcmpi(aud_lng_rgh(:,1),tbl_reg{iR}));
        if ~isempty(rgh_ind)
            tbl_tbl{iR,21} = aud_lng_rgh{rgh_ind,2};
        else
            tbl_tbl{iR,21} = 0;
        end
        
        rgh_ind = find(strcmpi(aud_con_rgh(:,1),tbl_reg{iR}));
        if ~isempty(rgh_ind)
            tbl_tbl{iR,22} = aud_con_rgh{rgh_ind,2};
        else
            tbl_tbl{iR,22} = 0;
        end
        
        rgh_ind = find(strcmpi(aud_phn_rgh(:,1),tbl_reg{iR}));
        if ~isempty(rgh_ind)
            tbl_tbl{iR,23} = aud_phn_rgh{rgh_ind,2};
        else
            tbl_tbl{iR,23} = 0;
        end
        
        rgh_ind = find(strcmpi(inc_rgh(:,1),tbl_reg{iR}));
        if ~isempty(rgh_ind)
            tbl_tbl{iR,25} = inc_rgh{rgh_ind,2};
        else
            tbl_tbl{iR,25} = 0;
        end
        
    else
        tbl_tbl{iR,15} = '-';
        
        tbl_tbl{iR,17} = '-';
        tbl_tbl{iR,18} = '-';
        tbl_tbl{iR,19} = '-';
        
        tbl_tbl{iR,21} = '-';
        tbl_tbl{iR,22} = '-';
        tbl_tbl{iR,23} = '-';
        
        tbl_tbl{iR,25} = '-';
    end
    
end

tbl_tbl = mmil_order_table2(tbl_tbl);

rmv_ind = [];
for iTB = 1:size(tbl_tbl,1)
    if ~isempty(tbl_tbl{iTB,1}) && ~ismember(tbl_tbl{iTB,1},inc_reg)
        rmv_ind = [rmv_ind iTB];
    end
end
tbl_tbl(rmv_ind,:) = [];

tbl_rmv_ind = cellfun(@(x) find(strcmpi(tbl_tbl(:,1),x)),{'Caudal Middle Frontal'},'uni',0);
tbl_tbl(tbl_rmv_ind{end}(2),:) = [];

%% If combining, then combine!
tbl_tbl_cmb = cell(size(row_cmb,1),21);

for iR = 1:numel(row_cmb)
    
    lft_ind = cellfun(@(x) find(strcmpi(tbl_tbl(:,1),x)),row_cmb{iR},'uni',0); lft_ind = [lft_ind{:}];
    
    if ~isempty(tbl_tbl(iR,1))
        
        tbl_tbl_cmb{iR,1} = row_nme{iR};
        
        for iC = 2:size(tbl_tbl,2)
            if ~isempty(tbl_tbl{1,iC})
                
                tbl_tbl_cmb{iR,iC} = sum([tbl_tbl{lft_ind,iC}]);
                
            end            
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmil_chk_dir('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/manuscript/tables')

cell2csv('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/manuscript/tables/table2.csv',[ovr_tbl_lbl(:,1:13)      ; tbl_lbl(:,1:13)      ; tbl_tbl(:,1:13)])
cell2csv('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/manuscript/tables/table3.csv',[ovr_tbl_lbl(:,[1 15:25]) ; tbl_lbl(:,[1 15:25]) ; tbl_tbl(:,[1 15:25])])

cell2csv('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/manuscript/tables/table2_cmb.csv',[ovr_tbl_lbl(:,1:13)      ; tbl_lbl(:,1:13)      ; tbl_tbl_cmb(:,1:13)])
cell2csv('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/manuscript/tables/table3_cmb.csv',[ovr_tbl_lbl(:,[1 15:25]) ; tbl_lbl(:,[1 15:25]) ; tbl_tbl_cmb(:,[1 15:25])])


end