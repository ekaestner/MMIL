paths = genpath('/home/ekaestne/gitrep/MMIL');
paths = strsplit(paths, pathsep);
filteredPaths = paths(~contains(paths, 'matlablt20')); % Example: Exclude 'matlablt20'
addpath(strjoin(filteredPaths, pathsep));

%%

bid_dir = '/space/mcdonald-syn01/1/BIDS/CAPES';

cov_dir = '/space/mcdonald-syn01/1/data/ai/original_covariates/';
key_fle = 'key_fle.csv';

prj_nme = 'CAPES';
cov_fle = 'CAPES_DATA_2025-04-28_1924_pipes.csv';
nme_fle = 'CAPES_BIDS_ejk.txt';

msc_dir = '/space/mcdonald-syn01/1/onboarding/locipipe_onboarding/MUSC_bids';
emo_dir = '/space/mcdonald-syn01/1/onboarding/locipipe_onboarding/EMORY_bids';

%%
fcfg = [];
fcfg.dta_loc = [cov_dir '/' key_fle];
[ key_dta, key_sbj, key_col ] = ejk_dta_frm(fcfg);

fcfg = [];
fcfg.dta_loc = [cov_dir '/' prj_nme '/' cov_fle];
[ cov_dta, cov_sbj, cov_col ] = ejk_dta_frm(fcfg);

fcfg = [];
fcfg.dta_loc = [cov_dir '/' prj_nme '/' nme_fle];
[ nme_dta, nme_sbj, nme_col ] = ejk_dta_frm(fcfg);

%% Pull relevant Redcap rows
rel_row = string_find(cov_dta(:,strcmpi(cov_col,'redcap_event_name')),'baseline');
cov_sbj = cov_sbj(rel_row,:);
cov_dta = cov_dta(rel_row,:);

%% Get CAPES subjects
bid_sbj = dir(bid_dir);
bid_sbj = {bid_sbj(:).name};
bid_sbj = bid_sbj(string_find(bid_sbj,{'sub-'}));

%% Fix Names normally
tot_rnm_fnd = 0;
cov_sbj_fix = cov_sbj;
for iS = 1:numel(nme_sbj)
    row_rpl = strcmpi( cov_sbj_fix, nme_sbj{iS});
    if any(row_rpl)
        cov_sbj_fix{row_rpl} = nme_dta{iS,1};
        tot_rnm_fnd = tot_rnm_fnd + 1;
    end
end

%% Fix BID's naming conventions
% cov_sbj_fix = cov_sbj;
% cov_sbj_fix = strrep(cov_sbj_fix,'_L_','_P_');
% cov_sbj_fix = strrep(cov_sbj_fix,'_R_','_P_');
% cov_sbj_fix = strrep(cov_sbj_fix,'_','');
% cov_sbj_fix = strrep(cov_sbj_fix,'EMU','EMO');
% cov_sbj_fix = strrep(cov_sbj_fix,'emu','EMO');
% cov_sbj_fix = strrep(cov_sbj_fix,'NTW','NOR');
cov_sbj_fix = strrep(cov_sbj_fix,'sub-','');

%% Find matching rows in main file
fnd_sbj = [];
mss_red_sbj = [ cov_sbj cov_sbj_fix ];
mss_red_rmv = zeros(size(mss_red_sbj,1),1);
mss_bid_sbj = [];
dbl_sbj = cell(0,5);
for iS = 1:numel(bid_sbj)

    row_sch = find(strcmpi( cov_sbj_fix, bid_sbj{iS}(5:end) ));

    if ~isempty(row_sch)
        if numel(row_sch)==1
            fnd_sbj = [ fnd_sbj ; cov_sbj(row_sch) cov_sbj_fix(row_sch) bid_sbj(iS) ];
            mss_red_rmv(row_sch,1) = 1;
        else
            dbl_sbj = [ dbl_sbj ; cov_sbj(row_sch(1)) cov_sbj_fix(row_sch(1)) cov_sbj(row_sch(2)) cov_sbj_fix(row_sch(2)) bid_sbj(iS) ];
        end
    else
        mss_bid_sbj = [ mss_bid_sbj ; bid_sbj(iS) ];
    end
end
mss_red_sbj(logical(mss_red_rmv),:) = [];

cell2csv([ cov_dir '/' prj_nme '/' 'output' '/' 'found_subjects.csv'],fnd_sbj)
cell2csv([ cov_dir '/' prj_nme '/' 'output' '/' 'missing_redcap_subjects.csv'],mss_red_sbj)
cell2csv([ cov_dir '/' prj_nme '/' 'output' '/' 'missing_BIDS_subjects.csv'],mss_bid_sbj)
cell2csv([ cov_dir '/' prj_nme '/' 'output' '/' 'double_subjects.csv'],dbl_sbj)

cell2csv([ cov_dir '/' prj_nme '/' 'output' '/' 'bids_subjects.csv'],bid_sbj)
cell2csv([ cov_dir '/' prj_nme '/' 'output' '/' 'redcap_subjets.csv'],cov_sbj)

%% Doublecheck Redcap missing
sbj_nme_out = mss_red_sbj;
sbj_nme_out(string_find(sbj_nme_out(:,2),'Missing'),:) = [];
sbj_nme_out(string_find(sbj_nme_out(:,2),'double'),:) = [];
sbj_nme_out(string_find(sbj_nme_out(:,2),'unusable'),:) = [];

cov_sbj_ind = find(ismember( cov_sbj, sbj_nme_out(:,1)));

srg_col = [ string_find(cov_col,'ep_surg') ...
            string_find(cov_col,'retro_surg') ...
            string_find(cov_col,'retro_cand') ...
            string_find(cov_col,'sx_procedure') ...
            string_find(cov_col,'sx_procedure_oth') ];

srg_col_nme = cov_col(srg_col);
srg_dta_out = cov_dta(cov_sbj_ind,srg_col);

cell2csv([ cov_dir '/' prj_nme '/' 'output' '/' 'missing_redcap_subjects_surgical_investigation.csv'],[ {'RedcapID'} {'BIDS_ID'} srg_col_nme ; sbj_nme_out srg_dta_out ]);

%% Check Found Surgery

cov_sbj_ind = find(ismember( cov_sbj, fnd_sbj(:,1)));
srg_dta_fnd = cov_dta(cov_sbj_ind,srg_col);
sum(cell2mat(srg_dta_fnd(:,4:8)))


%% Assess BIDS missing
ste_nme = cellfun(@(x) x(5:7),mss_bid_sbj,'uni',0);
out_rpt = tabulate(ste_nme);

mss_bid_sbj_out = [ mss_bid_sbj num2cell(zeros(numel(mss_bid_sbj),1))];

% Check MUSC
msc_sbj = dir(msc_dir);
msc_sbj = {msc_sbj(:).name};
msc_sbj = msc_sbj(string_find(msc_sbj,{'sub-'}));
msc_sbj = strrep(msc_sbj,'MUSPR','MUSP');
msc_sbj = strrep(msc_sbj,'MUSPL','MUSP');

for iS = 1:numel(msc_sbj)
    sbj_ind = strcmpi(mss_bid_sbj_out(:,1),msc_sbj{iS});
    if any(sbj_ind)
        mss_bid_sbj_out(sbj_ind,2) = {1};
    end
end

% Check Emory
emo_sbj = dir(emo_dir);
emo_sbj = {emo_sbj(:).name};
emo_sbj = emo_sbj(string_find(emo_sbj,{'sub-'}));
emo_sbj = strrep(emo_sbj,'Emory','');

for iS = 1:numel(emo_sbj)
    sbj_ind = strcmpi(mss_bid_sbj_out(:,1),emo_sbj{iS});
    if any(sbj_ind)
        mss_bid_sbj_out(sbj_ind,2) = {1};
    end
end

% Final Check
for iR = 1:size(out_rpt,1)
    out_rpt{iR,3} = sum(cell2mat(mss_bid_sbj_out(strcmpi(ste_nme,out_rpt{iR,1}),2))==1);
    out_rpt{iR,4} = sum(cell2mat(mss_bid_sbj_out(strcmpi(ste_nme,out_rpt{iR,1}),2))==0);
end

[ {'Site'} {'BIDS not in CAPES Redcap'} {'Confirmed in ENIGMA'} {'NOT in ENIGMA'} ; out_rpt]

mss_bid_sbj_out(strcmpi(ste_nme,out_rpt{2,1}) & cell2mat(mss_bid_sbj_out(:,2))==0,1)

%% New renaming
% 643 CAPES subjects in Redcap
% 934 CAPES subject scans on server
% 598 'Renaming Rows' alledgedly
%
% 508 subjects matching
% 434 Missing BIDS subjects
% 135 Missing Redcap subjects
%
% 22 missing from site
% 5 unusable
% 108 missing

numel(cov_sbj)
numel(bid_sbj)
numel(nme_sbj)

size(fnd_sbj,1)
size(mss_bid_sbj,1)
size(mss_red_sbj,1)

numel(string_find(mss_red_sbj,'Missing_get_from_site'))
numel(string_find(mss_red_sbj,'Remove_unusable'))
size(mss_red_sbj,1) - numel(string_find(mss_red_sbj(:,2),'_'))

%% Without Renaming
% 643 CAPES subjects in Redcap
% 934 CAPES subject scans on server
% 342 'Renaming Rows' alledgedly
%
% 267 subjects matching
% 667 Missing BIDS subjects
% 376 Missing Redcap subjects

numel(cov_sbj)
numel(bid_sbj)
numel(nme_sbj)

size(fnd_sbj,1)
size(mss_bid_sbj,1)
size(mss_red_sbj,1)

pit_mtc = fnd_sbj(string_find(fnd_sbj(:,3),'PIT'),3);
pit_mss_bid = mss_bid_sbj(string_find(mss_bid_sbj,'PIT'));
pit_mss_red = mss_red_sbj(string_find(mss_red_sbj(:,2),'PIT'),2);
pit_org_nme = nme_sbj(string_find(nme_sbj,'PIT'));
pit_new_nme = nme_dta(string_find(nme_sbj,'PIT'));

%% With Renaming
% 643 CAPES subjects in Redcap
% 934 CAPES subject scans on server
% 342 'Renaming Rows' alledgedly
%
% 372 subjects matching
% 524 Missing BIDS subjects
% 271 Missing Redcap subjects
% 38 subjects with Double match

numel(cov_sbj)
numel(bid_sbj)
numel(nme_sbj)

size(fnd_sbj,1)
size(mss_bid_sbj,1)
size(mss_red_sbj,1)
size(dbl_sbj,1)

numel(string_find(fnd_sbj(:,3),'P0'))
numel(string_find(fnd_sbj(:,3),'C0'))

numel(string_find(mss_red_sbj(:,2),'P0'))
numel(string_find(mss_red_sbj(:,2),'C0'))

pit_mtc = fnd_sbj(string_find(fnd_sbj(:,3),'PIT'),3);
pit_mss_bid = mss_bid_sbj(string_find(mss_bid_sbj,'PIT'));
pit_mss_red = mss_red_sbj(string_find(mss_red_sbj(:,2),'PIT'),2);
pit_dbl_nme = dbl_sbj(string_find(dbl_sbj(:,1),'PIT'),:);

pen_mtc = fnd_sbj(string_find(fnd_sbj(:,3),'PEN'),3);
pen_mss_bid = mss_bid_sbj(string_find(mss_bid_sbj,'PEN'));
pen_mss_red = mss_red_sbj(string_find(mss_red_sbj(:,2),'PEN'),2);
pen_dbl_nme = dbl_sbj(string_find(dbl_sbj(:,1),'PEN'),:);

%%
% Found subjects
fnd_sbj_rel = fnd_sbj(~isnan(fnd_sbj));
numel(fnd_sbj_rel)
cov_sbj_fix(fnd_sbj_rel)

% Missing BIDS subjects
mss_sbj_rel = mss_sbj(~isnan(mss_sbj));
numel(mss_sbj_rel)
bid_sbj(mss_sbj_rel)

% 
mss_org_sbj = cov_sbj_fix;
mss_org_sbj(org_sbj) = [];

[ nme_dta(strcmpi(nme_sbj,'EMU_P_0003')) nme_dta(strcmpi(nme_sbj,'EMU_R_0003'))]

[ ]

ttt = tabulate(cov_sbj_fix);
[~,iR] = sort(cell2mat(ttt(:,2)));
ttt(iR,:)

[  nme_sbj nme_dta  ]


strrep(bid_sbj(string_find(bid_sbj,'PIT')),'sub-','')'

sort(cov_sbj_fix(string_find(cov_sbj_fix,'PIT')))

nme_sbj(string_find(nme_sbj,'PIT'))
nme_dta(string_find(nme_sbj,'PIT'))

find(strcmpi(bid_sbj,'sub-PITP0001'))


find(strcmpi(bid_sbj,'sub-EMOP0033'))
find(strcmpi(cov_sbj_fix,'EMOP0033'))
cov_sbj(strcmpi(cov_sbj_fix,'EMOP0033'))
[ nme_sbj(strcmpi(nme_sbj,'EMU_P_0019')) nme_dta(strcmpi(nme_sbj,'EMU_P_0019')) ]


