clear; clc;

prj_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Collaborations/Adam/epd_t2h/';
dta_dir = [ prj_dir '/' 'Data' '/'];
plt_dir = [ prj_dir '/' 'Output' '/'];

dta_fle = 'ashs_combined_group_only.csv';

%%
grp_nme = {  'Y-HC' 'Y-TLE' 'O-HC' 'O-Epilepsy' 'MCI-C'};
grp_pos = [  1      2       4      5            6];
grp_col = { rgb('light grey') rgb('light orange') rgb('dark grey') rgb('dark orange') rgb('dark purple') };

%%
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' dta_fle ];
[ hip_sub_dta, hip_sub_sbj, hip_sub_col ] = ejk_dta_frm(fcfg);

% unique(hip_sub_dta(:,strcmpi(hip_sub_col,'Group')));

hip_sub_dta( strcmpi(hip_sub_dta(:,strcmpi(hip_sub_col,'Group')),'HC') & cell2mat(hip_sub_dta(:,strcmpi(hip_sub_col,'Age')))>=55, strcmpi(hip_sub_col,'Group')) = {'O-HC'};
hip_sub_dta( strcmpi(hip_sub_dta(:,strcmpi(hip_sub_col,'Group')),'HC') & cell2mat(hip_sub_dta(:,strcmpi(hip_sub_col,'Age')))<55, strcmpi(hip_sub_col,'Group')) = {'Y-HC'};

plt_fld = hip_sub_col;

%%
for iF = 1:numel(plt_fld)
    
    if isnumeric(hip_sub_dta{1,iF})
        
        fcfg = [];
        
        for iG = 1:numel(grp_nme)
            fcfg.xdt{iG}         = grp_pos(iG);
            fcfg.ydt{iG}         = cell2mat(hip_sub_dta( strcmpi(hip_sub_dta(:,strcmpi(hip_sub_col,'Group')),grp_nme{iG}), strcmpi(hip_sub_col,plt_fld{iF}) ));
            fcfg.fce_col{iG}     = grp_col{iG};
            fcfg.box_plt_col{iG} = grp_col{iG};
        end
        
        fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.xdt)) ];
        fcfg.box_plt = ones(1,numel(fcfg.xdt));
        
        fcfg.xlm = [ 0.5 max(cell2mat(fcfg.xdt))+0.5 ];
        fcfg.xlb = grp_nme;
        
        fcfg.ylb = plt_fld(iF);
        
        fcfg.mkr_sze = repmat(20,1,numel(fcfg.xdt));
        fcfg.aph_val = 0.45;
        
        fcfg.out_dir = plt_dir;
        fcfg.out_nme = plt_fld{iF};
        
        ejk_scatter(fcfg)
        
    end
end




