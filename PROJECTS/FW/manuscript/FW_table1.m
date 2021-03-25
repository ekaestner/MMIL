function FW_table1

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
         
%% Get Data for Selective Electrodes  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical';
nme = 'pap_rsp_600';

sig_lef = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_lhs_table_plot']); sig_lbl = sig_lef(1,:);
sig_lef = sig_lef(2:end-1,1:6);
sig_lef(1:end,1) = cellfun(@(x) x(5:end),sig_lef(1:end,1),'uni',0);

sig_rgh = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_rhs_table_plot']);
sig_rgh = sig_rgh(2:end-1,1:6);
sig_rgh(1:end,1) = cellfun(@(x) x(5:end),sig_rgh(1:end,1),'uni',0);

rmv_ind = [];
for iR = 1:size(sig_lef,1)
    if sig_lef{iR,5}<2
        rmv_ind = [rmv_ind iR];
    end
end
% sig_lef(rmv_ind,:) = [];
sig_lef = mmil_order_table(sig_lef);
sig_lef = mmil_order_table(sig_lef);

rmv_ind = [];
for iR = 1:size(sig_rgh,1)
    if sig_rgh{iR,5}<2
        rmv_ind = [rmv_ind iR];
    end
end
% sig_rgh(rmv_ind,:) = [];
sig_rgh = mmil_order_table(sig_rgh);
sig_rgh = mmil_order_table(sig_rgh);

% Make Table For Selective Electrodes
tbl_reg = unique([sig_lef(:,1) ; sig_rgh(:,1)]);
% tbl_reg = mmil_order_table(tbl_reg);
tbl_reg(strcmpi(tbl_reg,'no_location')) = [];
tbl_reg(strcmpi(tbl_reg,'unknown')) = [];

ovr_tbl_lbl = [repmat({''},1,2) {'Left Hemisphere'} repmat({''},1,2) {'Right Hemisphere'} repmat({''},1,1)];
tbl_lbl     = {'Region' '' 'Language-Selective Electrodes' 'Subjects' '' 'Language-Selective Electrodes' 'Subjects'};
tbl_tbl     = cell(size(tbl_reg,1),7);

for iR = 1:numel(tbl_reg)

    lft_ind = find(strcmpi(sig_lef(:,1),tbl_reg{iR}));
    rgh_ind = find(strcmpi(sig_rgh(:,1),tbl_reg{iR}));
    
    tbl_tbl{iR,1} = tbl_reg{iR};
    
    if ~isempty(lft_ind) > 0
        tbl_tbl{iR,3} = [num2str(round((sig_lef{lft_ind,2}/sig_lef{lft_ind,3})*100)) '% (' num2str(sig_lef{lft_ind,2}) ' / ' num2str(sig_lef{lft_ind,3}) ')'];;
        tbl_tbl{iR,4} = [num2str(sig_lef{lft_ind,5}) ' / ' num2str(sig_lef{lft_ind,6})];
    else
        tbl_tbl{iR,3} = '-';
        tbl_tbl{iR,4} = '-';
    end
    
    if ~isempty(rgh_ind) > 0
        tbl_tbl{iR,6}  = [num2str(round((sig_rgh{rgh_ind,2}/sig_rgh{rgh_ind,3})*100)) '% (' num2str(sig_rgh{rgh_ind,2}) ' / ' num2str(sig_rgh{rgh_ind,3}) ')'];
        tbl_tbl{iR,7}  = [num2str(sig_rgh{rgh_ind,5}) ' / ' num2str(sig_rgh{rgh_ind,6})];
    else
        tbl_tbl{iR,6} = '-';
        tbl_tbl{iR,7} = '-';
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
loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical';
nme = 'pap_rsp_600';

sig_lef = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_lhs_table_plot']); sig_lbl = sig_lef(1,:);
sig_lef = sig_lef(2:end-1,1:6);
sig_lef(1:end,1) = cellfun(@(x) x(5:end),sig_lef(1:end,1),'uni',0);

sig_lef_nme = load([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_lhs_table_plot_subject_names.mat']);
sig_lef_nme.ttt_sav = sig_lef_nme.ttt_sav(2:end-1,1:6);
sig_lef_nme.ttt_sav(1:end,1) = cellfun(@(x) x(5:end),sig_lef_nme.ttt_sav(1:end,1),'uni',0);

sig_rgh = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_rhs_table_plot']);
sig_rgh = sig_rgh(2:end-1,1:6);
sig_rgh(1:end,1) = cellfun(@(x) x(5:end),sig_rgh(1:end,1),'uni',0);

sig_rgh_nme = load([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme '/' 'total' '/' 'total_' nme '_rhs_table_plot_subject_names.mat']);
sig_rgh_nme.ttt_sav = sig_rgh_nme.ttt_sav(2:end-1,1:6);
sig_rgh_nme.ttt_sav(1:end,1) = cellfun(@(x) x(5:end),sig_rgh_nme.ttt_sav(1:end,1),'uni',0);

rmv_ind = [];
for iR = 1:size(sig_lef,1)
    if sig_lef{iR,5}<2
        rmv_ind = [rmv_ind iR];
    end
end
% sig_lef(rmv_ind,:) = [];
sig_lef = mmil_order_table(sig_lef);
sig_lef = mmil_order_table(sig_lef);

rmv_ind = [];
for iR = 1:size(sig_rgh,1)
    if sig_rgh{iR,5}<2
        rmv_ind = [rmv_ind iR];
    end
end
% sig_rgh(rmv_ind,:) = [];
sig_rgh = mmil_order_table(sig_rgh);
sig_rgh = mmil_order_table(sig_rgh);

% Make Table For Selective Electrodes
tbl_reg = unique([sig_lef(:,1) ; sig_rgh(:,1)]);
% tbl_reg = mmil_order_table(tbl_reg);
tbl_reg(strcmpi(tbl_reg,'no_location')) = [];
tbl_reg(strcmpi(tbl_reg,'unknown')) = [];

ovr_tbl_lbl = [repmat({''},1,2) {'Left Hemisphere'} repmat({''},1,2) {'Right Hemisphere'} repmat({''},1,1)];
tbl_lbl     = {'Region' '' 'Language-Selective Electrodes' 'Subjects' '' 'Language-Selective Electrodes' 'Subjects'};
tbl_tbl_cmb     = cell(size(row_cmb,1),7);

for iR = 1:numel(row_cmb)
    if ~isempty(row_cmb{iR})
    
    lft_ind = cellfun(@(x) find(strcmpi(sig_lef(:,1),x)),row_cmb{iR});
        lft_ind_nme = cellfun(@(x) find(strcmpi(sig_lef_nme.ttt_sav(:,1),x)),row_cmb{iR});
    rgh_ind = cellfun(@(x) find(strcmpi(sig_rgh(:,1),x)),row_cmb{iR});
        rgh_ind_nme = cellfun(@(x) find(strcmpi(sig_rgh_nme.ttt_sav(:,1),x)),row_cmb{iR});
    
    tbl_tbl_cmb{iR,1} = row_nme{iR};
    
    if ~isempty(lft_ind) > 0
        tbl_tbl_cmb{iR,3} = [num2str(round(sum([sig_lef{lft_ind,2}]) / sum([sig_lef{lft_ind,3}]) * 100)) '% (' num2str(sum([sig_lef{lft_ind,2}])) ' / ' num2str(sum([sig_lef{lft_ind,3}])) ')'];
        rsp_nme = unique(cat(1,sig_lef_nme.ttt_sav{lft_ind_nme,5}));
        tot_nme = unique(cat(1,sig_lef_nme.ttt_sav{lft_ind_nme,6}));
        tbl_tbl_cmb{iR,4} = [num2str(numel(rsp_nme)) ' / ' num2str(numel(tot_nme))];
    else
        tbl_tbl_cmb{iR,3} = '-';
        tbl_tbl_cmb{iR,4} = '-';
    end
    
    if ~isempty(rgh_ind) > 0
        tbl_tbl_cmb{iR,6}  = [num2str(round(sum([sig_rgh{rgh_ind,2}]) / sum([sig_rgh{rgh_ind,3}]) * 100)) '% (' num2str(sum([sig_rgh{rgh_ind,2}])) ' / ' num2str(sum([sig_rgh{rgh_ind,3}])) ')'];
        rsp_nme = unique(cat(1,sig_rgh_nme.ttt_sav{rgh_ind_nme,5}));
        tot_nme = unique(cat(1,sig_rgh_nme.ttt_sav{rgh_ind_nme,6}));
        tbl_tbl_cmb{iR,7} = [num2str(numel(rsp_nme)) ' / ' num2str(numel(tot_nme))];
    else
        tbl_tbl_cmb{iR,6} = '-';
        tbl_tbl_cmb{iR,7} = '-';
    end
    
    end
end

%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmil_chk_dir('/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/tables')

cell2csv('/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/tables/table1.csv',    [ovr_tbl_lbl ; tbl_lbl ; cell(1,numel(tbl_lbl)) ; tbl_tbl])
cell2csv('/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/tables/table1_cmb.csv',[ovr_tbl_lbl ; tbl_lbl ; cell(1,numel(tbl_lbl)) ; tbl_tbl_cmb])

%%
loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical';
nme = {'pap_rsp_600' 'pap_wrd_600' 'pap_lex_600' 'pap_con_600'};

sig_lef = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{1} '/' 'total' '/' 'total_' nme{1} '_lhs_table_plot']); sig_lbl = sig_lef(1,:);
sig_lef = sig_lef(2:end-1,1:6);
sig_lef(1:end,1) = cellfun(@(x) x(5:end),sig_lef(1:end,1),'uni',0);

wrd_lef = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{2} '/' 'total' '/' 'total_' nme{2} '_lhs_table_plot']); sig_lbl = wrd_lef(1,:);
wrd_lef = wrd_lef(2:end-1,[1:6 12:16]);
wrd_lef(1:end,1) = cellfun(@(x) x(5:end),wrd_lef(1:end,1),'uni',0);

lex_lef = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{3} '/' 'total' '/' 'total_' nme{3} '_lhs_table_plot']); sig_lbl = lex_lef(1,:);
lex_lef = lex_lef(2:end-1,:);
lex_lef(1:end,1) = cellfun(@(x) x(5:end),lex_lef(1:end,1),'uni',0);

con_lef = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{4} '/' 'total' '/' 'total_' nme{4} '_lhs_table_plot']); sig_lbl = con_lef(1,:);
con_lef = con_lef(2:end-1,1:6);
con_lef(1:end,1) = cellfun(@(x) x(5:end),con_lef(1:end,1),'uni',0);

sig_rgh = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{1} '/' 'total' '/' 'total_' nme{1} '_rhs_table_plot']);
sig_rgh = sig_rgh(2:end-1,1:6);
sig_rgh(1:end,1) = cellfun(@(x) x(5:end),sig_rgh(1:end,1),'uni',0);

wrd_rgh = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{2} '/' 'total' '/' 'total_' nme{2} '_rhs_table_plot']); sig_lbl = wrd_rgh(1,:);
wrd_rgh = wrd_rgh(2:end-1,[1:6 12:16]);
wrd_rgh(1:end,1) = cellfun(@(x) x(5:end),wrd_rgh(1:end,1),'uni',0);

lex_rgh = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{3} '/' 'total' '/' 'total_' nme{3} '_rhs_table_plot']); sig_lbl = lex_rgh(1,:);
lex_rgh = lex_rgh(2:end-1,:);
lex_rgh(1:end,1) = cellfun(@(x) x(5:end),lex_rgh(1:end,1),'uni',0);

con_rgh = mmil_readtext([loc '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' nme{4} '/' 'total' '/' 'total_' nme{4} '_rhs_table_plot']); sig_lbl = con_rgh(1,:);
con_rgh = con_rgh(2:end-1,1:6);
con_rgh(1:end,1) = cellfun(@(x) x(5:end),con_rgh(1:end,1),'uni',0);

rmv_ind = [];
for iR = 1:size(sig_lef,1)
    if sig_lef{iR,5}<2
        rmv_ind = [rmv_ind iR];
    end
end

% sig_lef(rmv_ind,:) = [];
sig_lef = mmil_order_table(sig_lef);
sig_lef = mmil_order_table(sig_lef);

% wrd_lef(rmv_ind,:) = [];
wrd_lef = mmil_order_table(wrd_lef);
wrd_lef = mmil_order_table(wrd_lef);

% lex_lef(rmv_ind,:) = [];
lex_lef = mmil_order_table(lex_lef);
lex_lef = mmil_order_table(lex_lef);

% con_lef(rmv_ind,:) = [];
con_lef = mmil_order_table(con_lef);
con_lef = mmil_order_table(con_lef);

rmv_ind = [];
for iR = 1:size(sig_rgh,1)
    if sig_rgh{iR,5}<2
        rmv_ind = [rmv_ind iR];
    end
end

% sig_rgh(rmv_ind,:) = [];
sig_rgh = mmil_order_table(sig_rgh);
sig_rgh = mmil_order_table(sig_rgh);

% wrd_rgh(rmv_ind,:) = [];
wrd_rgh = mmil_order_table(wrd_rgh);
wrd_rgh = mmil_order_table(wrd_rgh);

% lex_rgh(rmv_ind,:) = [];
lex_rgh = mmil_order_table(lex_rgh);
lex_rgh = mmil_order_table(lex_rgh);

% con_rgh(rmv_ind,:) = [];
con_rgh = mmil_order_table(con_rgh);
con_rgh = mmil_order_table(con_rgh);

% Make Table For Selective Electrodes
tbl_reg = unique([sig_lef(:,1) ; sig_rgh(:,1)]);
% tbl_reg = mmil_order_table(tbl_reg);
tbl_reg(strcmpi(tbl_reg,'no_location')) = [];
tbl_reg(strcmpi(tbl_reg,'unknown')) = [];

ovr_tbl_lbl = [repmat({''},1,7) {'Left Hemisphere'} repmat({''},1,9) {'Right Hemisphere'} repmat({''},1,3)];
tbl_lbl     = {'Region' '' 'Selective Electrodes' 'Orthographic-Selective Electrodes' 'Word-Selective Electrodes' '' 'Bigram-Frequency' 'Word-Frequency' 'Repetition' '' 'False-Font' '' 'Selective Electrodes' 'Orthographic-Selective Electrodes' 'Word-Selective Electrodes' '' 'Bigram-Frequency' 'Word-Frequency' 'Repetition' '' 'False-Font'};
tbl_tbl     = cell(size(tbl_reg,1),21);

for iR = 1:numel(tbl_reg)

    lft_ind = find(strcmpi(sig_lef(:,1),tbl_reg{iR}));
    rgh_ind = find(strcmpi(sig_rgh(:,1),tbl_reg{iR}));
    
    tbl_tbl{iR,1} = tbl_reg{iR};
    
    if ~isempty(lft_ind) > 0
        tbl_tbl{iR,3} = sig_lef{lft_ind,2};
        tbl_tbl{iR,4} = wrd_lef{lft_ind,2};
        tbl_tbl{iR,5} = wrd_lef{lft_ind,7};
        
        tbl_tbl{iR,7} = lex_lef{lft_ind,7};
        tbl_tbl{iR,8} = lex_lef{lft_ind,12};
        tbl_tbl{iR,9} = lex_lef{lft_ind,2};
        
        tbl_tbl{iR,11} = con_lef{lft_ind,2};
    else
        tbl_tbl{iR,3} = '-';
        tbl_tbl{iR,4} = '-';
        tbl_tbl{iR,5} = '-';
        
        tbl_tbl{iR,7} = '-';
        tbl_tbl{iR,8} = '-';
        tbl_tbl{iR,9} = '-';
        
        tbl_tbl{iR,11} = '-';
    end
    
    if ~isempty(rgh_ind) > 0
        tbl_tbl{iR,13} = sig_rgh{rgh_ind,2};
        tbl_tbl{iR,14} = wrd_rgh{rgh_ind,2};
        tbl_tbl{iR,15} = wrd_rgh{rgh_ind,7};
        
        tbl_tbl{iR,17} = lex_rgh{rgh_ind,7};
        tbl_tbl{iR,18} = lex_rgh{rgh_ind,12};
        tbl_tbl{iR,19} = lex_rgh{rgh_ind,2};
        
        tbl_tbl{iR,21} = con_rgh{rgh_ind,2};
    else
        tbl_tbl{iR,13} = '-';
        tbl_tbl{iR,14} = '-';
        tbl_tbl{iR,15} = '-';
        
        tbl_tbl{iR,17} = '-';
        tbl_tbl{iR,18} = '-';
        tbl_tbl{iR,19} = '-';
        
        tbl_tbl{iR,21} = '-';
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

%%

mmil_chk_dir('/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/tables')

cell2csv('/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/tables/table2.csv',[ovr_tbl_lbl(:,1:11) ; tbl_lbl(:,1:11) ; tbl_tbl(:,1:11)])
cell2csv('/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/tables/table3.csv',[ovr_tbl_lbl(:,[1 12:21]) ; tbl_lbl(:,[1 12:21]) ; tbl_tbl(:,[1 12:21])])

cell2csv('/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/tables/table2_cmb.csv',[ovr_tbl_lbl(:,1:11) ; tbl_lbl(:,1:11) ; tbl_tbl_cmb(:,1:11)])
cell2csv('/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/tables/table3_cmb.csv',[ovr_tbl_lbl(:,[1 12:21]) ; tbl_lbl(:,[1 12:21]) ; tbl_tbl_cmb(:,[1 12:21])])

end