clear; clc;

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
loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';

nme = 'pap_vis_act';
sig_lef{1} = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_lhs_table_plot']); sig_lbl = sig_lef{1}(1,:);
sig_lef{1} = sig_lef{1}(2:end-1,1:6);
sig_lef{1}(:,1) = cellfun(@(x) x(5:end),sig_lef{1}(:,1),'uni',0);

nme = 'pap_aud_act';
sig_lef{2} = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_lhs_table_plot']); sig_lbl = sig_lef{2}(1,:);
sig_lef{2} = sig_lef{2}(2:end-1,1:6);
sig_lef{2}(:,1) = cellfun(@(x) x(5:end),sig_lef{2}(:,1),'uni',0);

nme = 'pap_bim_act';
sig_lef{3} = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_lhs_table_plot']); sig_lbl = sig_lef{3}(1,:);
sig_lef{3} = sig_lef{3}(2:end-1,1:6);
sig_lef{3}(:,1) = cellfun(@(x) x(5:end),sig_lef{3}(:,1),'uni',0);

% Make Table For Selective Electrodes
tot_reg = unique([sig_lef{1}(:,1) ; sig_lef{2}(:,1) ; sig_lef{3}(:,1)]);
tbl_reg = mmil_order_table(tot_reg); tbl_reg = mmil_order_table(tbl_reg);
sig_lef{1} = mmil_order_table(sig_lef{1}); sig_lef{1} = mmil_order_table(sig_lef{1});
sig_lef{2} = mmil_order_table(sig_lef{2}); sig_lef{2} = mmil_order_table(sig_lef{2});
sig_lef{3} = mmil_order_table(sig_lef{3}); sig_lef{3} = mmil_order_table(sig_lef{3});

tbl_reg(strcmpi(tbl_reg,'no_location')) = [];
tbl_reg(strcmpi(tbl_reg,'unknown')) = [];

ovr_tbl_lbl = [repmat({''},1,2)  repmat({''},1,3)                               {'Left Hemisphere'} repmat({''},1,2)            repmat({''},1,2) ];
tbl_lbl     = {'Region' ''       'Visual Responsive Electrodes' 'Subjects' ''   'Auditory Responsive Electrodes' 'Subjects' ''  'Bi-Modal Responsive Electrodes' 'Subjects'};
tbl_tbl     = cell(size(tbl_reg,1),10);

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
    
    lft_ind = find(strcmpi(sig_lef{3}(:,1),tbl_reg{iR}));
    if ~isempty(lft_ind) > 0
        tbl_tbl{iR,9} = [num2str(round((sig_lef{3}{lft_ind,2}/sig_lef{3}{lft_ind,3})*100)) '% (' [num2str(sig_lef{3}{lft_ind,2}) ' / ' num2str(sig_lef{3}{lft_ind,3})] ')'];
        tbl_tbl{iR,10} = [num2str(sig_lef{3}{lft_ind,5}) ' / ' num2str(sig_lef{3}{lft_ind,6})];
    else
        tbl_tbl{iR,9} = '-';
        tbl_tbl{iR,10} = '-';
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

% If combining, then combine! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';

nme = 'pap_vis_act';
sig_lef{1} = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_lhs_table_plot']); sig_lbl = sig_lef{1}(1,:);
sig_lef{1} = sig_lef{1}(2:end-1,1:6);
sig_lef{1}(:,1) = cellfun(@(x) x(5:end),sig_lef{1}(:,1),'uni',0);

sig_lef_nme = load([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_lhs_table_plot_subject_names.mat']);
sig_lef_nme.ttt_sav = sig_lef_nme.ttt_sav(2:end-1,1:6);
sig_lef_nme.ttt_sav(1:end,1) = cellfun(@(x) x(5:end),sig_lef_nme.ttt_sav(1:end,1),'uni',0);
sig_lef_nme_hld{1} = sig_lef_nme.ttt_sav;

nme = 'pap_aud_act';
sig_lef{2} = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_lhs_table_plot']); sig_lbl = sig_lef{2}(1,:);
sig_lef{2} = sig_lef{2}(2:end-1,1:6);
sig_lef{2}(:,1) = cellfun(@(x) x(5:end),sig_lef{2}(:,1),'uni',0);

sig_lef_nme = load([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_lhs_table_plot_subject_names.mat']);
sig_lef_nme.ttt_sav = sig_lef_nme.ttt_sav(2:end-1,1:6);
sig_lef_nme.ttt_sav(1:end,1) = cellfun(@(x) x(5:end),sig_lef_nme.ttt_sav(1:end,1),'uni',0);
sig_lef_nme_hld{2} = sig_lef_nme.ttt_sav;

nme = 'pap_bim_act';
sig_lef{3} = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_lhs_table_plot']); sig_lbl = sig_lef{3}(1,:);
sig_lef{3} = sig_lef{3}(2:end-1,1:6);
sig_lef{3}(:,1) = cellfun(@(x) x(5:end),sig_lef{3}(:,1),'uni',0);

sig_lef_nme = load([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_lhs_table_plot_subject_names.mat']);
sig_lef_nme.ttt_sav = sig_lef_nme.ttt_sav(2:end-1,1:6);
sig_lef_nme.ttt_sav(1:end,1) = cellfun(@(x) x(5:end),sig_lef_nme.ttt_sav(1:end,1),'uni',0);
sig_lef_nme_hld{3} = sig_lef_nme.ttt_sav;

% Make Table For Selective Electrodes
tot_reg = unique([sig_lef{1}(:,1) ; sig_lef{2}(:,1) ; sig_lef{3}(:,1)]);
tbl_reg = mmil_order_table(tot_reg); tbl_reg = mmil_order_table(tbl_reg);
sig_lef{1} = mmil_order_table(sig_lef{1}); sig_lef{1} = mmil_order_table(sig_lef{1});
sig_lef{2} = mmil_order_table(sig_lef{2}); sig_lef{2} = mmil_order_table(sig_lef{2});
sig_lef{3} = mmil_order_table(sig_lef{3}); sig_lef{3} = mmil_order_table(sig_lef{3});

tbl_reg(strcmpi(tbl_reg,'no_location')) = [];
tbl_reg(strcmpi(tbl_reg,'unknown')) = [];

ovr_tbl_lbl = [repmat({''},1,2)  repmat({''},1,3)                               {'Left Hemisphere'} repmat({''},1,2)            repmat({''},1,2) ];
tbl_lbl     = {'Region' ''       'Visual Responsive Electrodes' 'Subjects' ''   'Auditory Responsive Electrodes' 'Subjects' ''  'Bi-Modal Responsive Electrodes' 'Subjects'};
tbl_tbl_cmb     = cell(size(row_cmb,1),10);

for iR = 1:numel(row_nme)
    
    tbl_tbl_cmb{iR,1} = row_nme{iR};
    
    lft_ind = cellfun(@(x) find(strcmpi(sig_lef{1}(:,1),x)),row_cmb{iR},'uni',0); lft_ind = [lft_ind{:}];
    lft_ind_nme = cellfun(@(x) find(strcmpi(sig_lef_nme_hld{1}(:,1),x)),row_cmb{iR});
    if ~isempty(lft_ind) > 0
        tbl_tbl_cmb{iR,3} = [num2str(round(( sum([sig_lef{1}{lft_ind,2}]) / sum([sig_lef{1}{lft_ind,3}]) )*100)) '% (' [num2str( sum([sig_lef{1}{lft_ind,2}]) ) ' / ' num2str( sum([sig_lef{1}{lft_ind,3}]) )] ')'];
        rsp_nme = unique(cat(1,sig_lef_nme_hld{1}{lft_ind_nme,5}));
        tot_nme = unique(cat(1,sig_lef_nme_hld{1}{lft_ind_nme,6}));
        tbl_tbl_cmb{iR,4} = [num2str(numel(rsp_nme)) ' / ' num2str(numel(tot_nme))];
    else
        tbl_tbl_cmb{iR,3} = '-';
        tbl_tbl_cmb{iR,4} = '-';
    end
    
    lft_ind = cellfun(@(x) find(strcmpi(sig_lef{3}(:,1),x)),row_cmb{iR},'uni',0); lft_ind = [lft_ind{:}];
    lft_ind_nme = cellfun(@(x) find(strcmpi(sig_lef_nme_hld{2}(:,1),x)),row_cmb{iR});
    if ~isempty(lft_ind) > 0
        tbl_tbl_cmb{iR,6} = [num2str(round(( sum([sig_lef{2}{lft_ind,2}]) / sum([sig_lef{2}{lft_ind,3}]) )*100)) '% (' [num2str( sum([sig_lef{2}{lft_ind,2}]) ) ' / ' num2str( sum([sig_lef{2}{lft_ind,3}]) )] ')'];
        rsp_nme = unique(cat(1,sig_lef_nme_hld{2}{lft_ind_nme,5}));
        tot_nme = unique(cat(1,sig_lef_nme_hld{2}{lft_ind_nme,6}));
        tbl_tbl_cmb{iR,7} = [num2str(numel(rsp_nme)) ' / ' num2str(numel(tot_nme))];
    else
        tbl_tbl_cmb{iR,6} = '-';
        tbl_tbl_cmb{iR,7} = '-';
    end
    
    lft_ind = cellfun(@(x) find(strcmpi(sig_lef{3}(:,1),x)),row_cmb{iR},'uni',0); lft_ind = [lft_ind{:}];
    lft_ind_nme = cellfun(@(x) find(strcmpi(sig_lef_nme_hld{3}(:,1),x)),row_cmb{iR});
    if ~isempty(lft_ind) > 0
        tbl_tbl_cmb{iR,9} = [num2str(round(( sum([sig_lef{3}{lft_ind,2}]) / sum([sig_lef{3}{lft_ind,3}]) )*100)) '% (' [num2str( sum([sig_lef{3}{lft_ind,2}]) ) ' / ' num2str( sum([sig_lef{3}{lft_ind,3}]) )] ')'];
        rsp_nme = unique(cat(1,sig_lef_nme_hld{3}{lft_ind_nme,5}));
        tot_nme = unique(cat(1,sig_lef_nme_hld{3}{lft_ind_nme,6}));
        tbl_tbl_cmb{iR,10} = [num2str(numel(rsp_nme)) ' / ' num2str(numel(tot_nme))];
    else
        tbl_tbl_cmb{iR,9} = '-';
        tbl_tbl_cmb{iR,10} = '-';
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmil_chk_dir('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/manuscript/tables')

cell2csv('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/manuscript/tables/table1.csv',[ovr_tbl_lbl ; tbl_lbl ; tbl_tbl])
cell2csv('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/manuscript/tables/table1_cmb.csv',[ovr_tbl_lbl ; tbl_lbl ; tbl_tbl_cmb])

%% Effects Table
loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';

% Check
nme = 'pap_vis_act';
sig_lef{1} = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_lhs_table_plot']); sig_lbl = sig_lef{1}(1,:);
sig_lef{1} = sig_lef{1}(2:end-1,1:6);
sig_lef{1}(:,1) = cellfun(@(x) x(5:end),sig_lef{1}(:,1),'uni',0);

nme = 'pap_aud_act';
sig_lef{2} = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_lhs_table_plot']); sig_lbl = sig_lef{2}(1,:);
sig_lef{2} = sig_lef{2}(2:end-1,1:6);
sig_lef{2}(:,1) = cellfun(@(x) x(5:end),sig_lef{2}(:,1),'uni',0);

nme = 'pap_bim_act';
sig_lef{3} = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_lhs_table_plot']); sig_lbl = sig_lef{3}(1,:);
sig_lef{3} = sig_lef{3}(2:end-1,1:6);
sig_lef{3}(:,1) = cellfun(@(x) x(5:end),sig_lef{3}(:,1),'uni',0);

% Load Left
nme = { 'pap_vis_rep' 'pap_vis_rep_nob' 'pap_aud_rep' 'pap_aud_rep_nob' 'pap_bim_rep' 'pap_bim_rep_nob'};

vis_rep_lef = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{1} '/' 'total' '/' 'total_' nme{1} '_lhs_table_plot']); sig_lbl = vis_rep_lef(1,:);
vis_rep_lef = vis_rep_lef(2:end-1,1:6);
vis_rep_lef(:,1) = cellfun(@(x) x(5:end),vis_rep_lef(:,1),'uni',0);

vis_frq_lef = mmil_readtext([loc '/' 'sig_chn' '/' 'lfp' '/' 'ecog' '/' 'split' '/' nme{2} '/' 'total' '/' 'total_' nme{2} '_lhs_table_plot']); sig_lbl = vis_frq_lef(1,:);
vis_frq_lef = vis_frq_lef(2:end-1,1:6);
vis_frq_lef(:,1) = cellfun(@(x) x(5:end),vis_frq_lef(:,1),'uni',0);

aud_rep_lef = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{3} '/' 'total' '/' 'total_' nme{3} '_lhs_table_plot']); sig_lbl = aud_rep_lef(1,:);
aud_rep_lef = aud_rep_lef(2:end-1,1:6);
aud_rep_lef(:,1) = cellfun(@(x) x(5:end),aud_rep_lef(:,1),'uni',0);

aud_frq_lef = mmil_readtext([loc '/' 'sig_chn' '/' 'lfp' '/' 'ecog' '/' 'split' '/' nme{4} '/' 'total' '/' 'total_' nme{4} '_lhs_table_plot']); sig_lbl = aud_frq_lef(1,:);
aud_frq_lef = aud_frq_lef(2:end-1,1:6);
aud_frq_lef(:,1) = cellfun(@(x) x(5:end),aud_frq_lef(:,1),'uni',0);

bim_rep_lef = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{5} '/' 'total' '/' 'total_' nme{5} '_lhs_table_plot']); sig_lbl = bim_rep_lef(1,:);
bim_rep_lef = bim_rep_lef(2:end-1,1:6);
bim_rep_lef(:,1) = cellfun(@(x) x(5:end),bim_rep_lef(:,1),'uni',0);

bim_frq_lef = mmil_readtext([loc '/' 'sig_chn' '/' 'lfp' '/' 'ecog' '/' 'split' '/' nme{6} '/' 'total' '/' 'total_' nme{6} '_lhs_table_plot']); sig_lbl = bim_frq_lef(1,:);
bim_frq_lef = bim_frq_lef(2:end-1,1:6);
bim_frq_lef(:,1) = cellfun(@(x) x(5:end),bim_frq_lef(:,1),'uni',0);

% Order
sig_lef{1} = mmil_order_table(sig_lef{1});
sig_lef{2} = mmil_order_table(sig_lef{2});
sig_lef{3} = mmil_order_table(sig_lef{3});

vis_rep_lef = mmil_order_table(vis_rep_lef);
vis_frq_lef = mmil_order_table(vis_frq_lef);

aud_rep_lef = mmil_order_table(aud_rep_lef);
aud_frq_lef = mmil_order_table(aud_frq_lef);

bim_rep_lef = mmil_order_table(bim_rep_lef);
bim_frq_lef = mmil_order_table(bim_frq_lef);

% Make Table For Selective Electrodes
tot_reg = unique([vis_rep_lef(:,1) ; vis_frq_lef(:,1) ; aud_rep_lef(:,1) ; aud_frq_lef(:,1) ; bim_rep_lef(:,1) ; bim_frq_lef(:,1)]);
if ~isempty(string_find(tot_reg,{' '})); tbl_reg = tot_reg; else tbl_reg = mmil_order_table(tot_reg); end
tbl_reg(strcmpi(tbl_reg,'no_location')) = [];
tbl_reg(strcmpi(tbl_reg,'unknown')) = [];

ovr_tbl_lbl = [repmat({''},1,4) {'Left Hemisphere'} repmat({''},1,5)];
tbl_lbl     = {'Region' '' 'Text Repetition' 'Text Repetition LFP' '' 'Voice Repetition' 'Voice Repetition LFP' '' 'Bimodal Repetition' 'Bimodal Repetition LFP'};
tbl_tbl     = cell(size(tbl_reg,1),10);

%
for iR = 1:numel(tbl_reg)
    
    tbl_tbl{iR,1} = tbl_reg{iR};
    
    % Left
    lft_ind_one = find(strcmpi(sig_lef{1}(:,1),tbl_reg{iR}));
    lft_ind_two = find(strcmpi(sig_lef{2}(:,1),tbl_reg{iR}));
    lft_ind_thr = find(strcmpi(sig_lef{3}(:,1),tbl_reg{iR}));
    
    if ~isempty(lft_ind_one) || ~isempty(lft_ind_two) || ~isempty(lft_ind_thr)
        
        lft_ind = find(strcmpi(vis_rep_lef(:,1),tbl_reg{iR}));
        if ~isempty(lft_ind)
            tbl_tbl{iR,3} = vis_rep_lef{lft_ind,2};
        else
            tbl_tbl{iR,3} = 0;
        end

        lft_ind = find(strcmpi(vis_frq_lef(:,1),tbl_reg{iR}));
        if ~isempty(lft_ind)
            tbl_tbl{iR,4} = vis_frq_lef{lft_ind,2};
        else
            tbl_tbl{iR,4} = 0;
        end
        
        lft_ind = find(strcmpi(aud_rep_lef(:,1),tbl_reg{iR}));
        if ~isempty(lft_ind)
            tbl_tbl{iR,6} = aud_rep_lef{lft_ind,2};
        else
            tbl_tbl{iR,6} = 0;
        end
        
        lft_ind = find(strcmpi(aud_frq_lef(:,1),tbl_reg{iR}));
        if ~isempty(lft_ind)
            tbl_tbl{iR,7} = aud_frq_lef{lft_ind,2};
        else
            tbl_tbl{iR,7} = 0;
        end
        
        
        lft_ind = find(strcmpi(bim_rep_lef(:,1),tbl_reg{iR}));
        if ~isempty(lft_ind)
            tbl_tbl{iR,9}  = bim_rep_lef{lft_ind,2};
        else
            tbl_tbl{iR,9} = 0;
        end
        
        
        lft_ind = find(strcmpi(bim_frq_lef(:,1),tbl_reg{iR}));
        if ~isempty(lft_ind)
            tbl_tbl{iR,10} = bim_frq_lef{lft_ind,2};
        else
            tbl_tbl{iR,10} = 0;
        end
        
    else
        tbl_tbl{iR,3} = '-';
        tbl_tbl{iR,4} = '-';
        
        tbl_tbl{iR,6} = '-';
        tbl_tbl{iR,7} = '-';
        
        tbl_tbl{iR,9}  = '-';
        tbl_tbl{iR,10} = '-';
    end
end

tbl_tbl = mmil_order_table(tbl_tbl);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tbl_tbl_cmb = cell(size(row_cmb,1),10);

for iR = 1:3 %numel(row_cmb)
    
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
mmil_chk_dir('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/manuscript/tables')

cell2csv('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/manuscript/tables/table2.csv',    [ovr_tbl_lbl(:,1:10)      ; tbl_lbl(:,1:10)      ; tbl_tbl(:,1:10)])
cell2csv('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/manuscript/tables/table2_cmb.csv',[ovr_tbl_lbl(:,1:10)      ; tbl_lbl(:,1:10)      ; tbl_tbl_cmb(:,1:10)])




























