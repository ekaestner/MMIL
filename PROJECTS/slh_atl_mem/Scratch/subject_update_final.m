clear; clc;

out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/slh_atl_mem/Analysis/subject_check';

%% Load pieces
% Original Emory
fcfg = [];
fcfg.dta_loc = [ out_dir '/' 'Emory_Check.csv'];
fcfg.dta_col = 2;
[ emy_org_dta, emy_org_sbj, emy_org_col] = ejk_dta_frm( fcfg );

% Meeting Notes Emory
fcfg = [];
fcfg.dta_loc = [ out_dir '/' 'RebeccaDanNotes.csv'];
fcfg.dta_col = 2;
[ emy_nte_dta, emy_nte_sbj, emy_nte_col] = ejk_dta_frm( fcfg );

% Final Sample
fcfg = [];
fcfg.dta_loc = [ out_dir '/' 'FinalSample.csv'];
fcfg.dta_col = 2;
[ fnl_smp_dta, fnl_smp_sbj, fnl_smp_col] = ejk_dta_frm( fcfg );

% Final Out
fcfg = [];
fcfg.dta_loc = [ out_dir '/' 'Final_Subjects_v1.csv'];
fcfg.dta_col = 2;
[ fnl_out_dta, fnl_out_sbj, fnl_out_col] = ejk_dta_frm( fcfg );

%% Check if any EPK neglected 
% from Dan's dataset %%%%%%%%%%%%%%
epk_sbj = fnl_smp_sbj(string_find(fnl_smp_sbj,'EPK'));

[ ~, ~, smp_ovr_lap] = intersect(fnl_out_sbj,epk_sbj);
mss_org_sbj          = epk_sbj(setxor(1:numel(epk_sbj),smp_ovr_lap));

[ ~, smp_ovr_lap, ~] = intersect(emy_org_sbj,mss_org_sbj);
[ emy_org_sbj(smp_ovr_lap) emy_org_dta(smp_ovr_lap, strcmpi(emy_org_col,'DanSurgery')) emy_org_dta(smp_ovr_lap, strcmpi(emy_org_col,'RebeccaSurgery'))];

%% Add Meeting Notes
add_col     = { 'RebeccaNotes'   'RebeccaSurgery'  'RebeccaReview' 'Notes from Meeting' 'Note from Email'};
add_col_nme = { 'RebeccaSurgery' 'SurgeryCategory' 'RebeccaReview' 'DanMeetingNotes'    'DanEmailNotes' };
add_tbl = cell(numel(fnl_out_sbj),numel(add_col));
for iS = 1:numel(fnl_out_sbj)
    sbj_ind = strcmpi(emy_nte_sbj,fnl_out_sbj{iS});
    for iC = 1:numel(add_col)
        col_ind = strcmpi(emy_nte_col,add_col{iC});
        add_tbl{iS,iC} = emy_nte_dta{sbj_ind,col_ind};
    end    
end

cell2csv([ out_dir '/' 'Final_Subjects_v2.csv'],[ 'sbj_nme' fnl_out_col add_col_nme ; fnl_out_sbj  fnl_out_dta add_tbl]);






