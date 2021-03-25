function mmil_add_subject(cfg)

% cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical';
% cfg.sbj_nme = 'NY007_SA_SZ';

if exist([cfg.clr_fld '/' 'subjects'],'file')
   sub_lst = mmil_readtext([cfg.clr_fld '/' 'subjects']);
else
   sub_lst = {cfg.sbj_nme};
   cell2csv([cfg.clr_fld '/' 'subjects'],sub_lst);
end

if ~any(strcmpi(cfg.sbj_nme,sub_lst)) 
   sub_lst = [sub_lst ; cfg.sbj_nme];
   cell2csv([cfg.clr_fld '/' 'subjects'],sub_lst);
end

sbj_num = find(strcmpi(sub_lst,cfg.sbj_nme));

if ~exist([cfg.clr_fld '/' 'sbj_inf' '/'],'dir'); mkdir([cfg.clr_fld '/' 'sbj_inf' '/']); end
if ~exist([cfg.clr_fld '/' 'sbj_inf' '/' cfg.sbj_nme],'file')

    bse_sbj_inf = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/bse_sbj_inf');
    for iR = 1:size(bse_sbj_inf,1)
        if ~isempty(bse_sbj_inf{iR}); sbj_inf_str{iR,1} = [bse_sbj_inf{iR} ' | []']; end
    end
    
    cell2csv([cfg.clr_fld '/' 'sbj_inf' '/' cfg.sbj_nme],sbj_inf_str);
    system(['gedit ' cfg.clr_fld '/' 'sbj_inf' '/' cfg.sbj_nme ' &']);
else
    clr_fle = mmil_readtext([cfg.clr_fld '/' 'sbj_inf' '/' cfg.sbj_nme],['`']);
    clr_fle = order_sbj_inf(clr_fle);
    cell2csv([cfg.clr_fld '/' 'sbj_inf' '/' cfg.sbj_nme],clr_fle);
    system(['gedit ' cfg.clr_fld '/' 'sbj_inf' '/' cfg.sbj_nme ' &']);
end

fprintf('\n\n%s Subject Number: %i\n\n',cfg.sbj_nme,sbj_num)

end
%% SUBFUNCTIONS



%% OLD
% sbj_inf_str{1,1} = 'bse_frq : []';
% sbj_inf_str{2,1} = 'cln_dir : []';
% sbj_inf_str{3,1} = 'cln_fld : []';
% sbj_inf_str{4,1} = 'cmn_nse : []';
% sbj_inf_str{5,1} = 'cmn_nme : {}';
% sbj_inf_str{6,1} = 'electrode_anatomical_location_depth : []';
% sbj_inf_str{7,1} = 'electrode_anatomical_location_ecog : []';
% sbj_inf_str{8,1} = 'electrode_average_location : []';
% sbj_inf_str{9,1} = 'electrode_location : []';
% sbj_inf_str{10,1} = 'end_dir : []';
% sbj_inf_str{11,1} = 'hemi : []';
% sbj_inf_str{12,1} = 'indir : []';
% sbj_inf_str{13,1} = 'pial : []';
% sbj_inf_str{14,1} = 'recon_location : []';
% sbj_inf_str{15,1} = 'rmv_chn : []';
% sbj_inf_str{16,1} = 'trialfun : []';
% sbj_inf_str{17,1} = 'tsk : []';
% sbj_inf_str{18,1} = 'ignore : []';

%% OLDER
% if ~exist([cfg.clr_fld '/' 'indir' '/'],'dir'); mkdir([cfg.clr_fld '/' 'indir' '/']); end
% if ~exist([cfg.clr_fld '/' 'cln_dir' '/'],'dir'); mkdir([cfg.clr_fld '/' 'cln_dir' '/']); end
% if ~exist([cfg.clr_fld '/' 'end_dir' '/'],'dir'); mkdir([cfg.clr_fld '/' 'end_dir' '/']); end
% if ~exist([cfg.clr_fld '/' 'tsk' '/'],'dir'); mkdir([cfg.clr_fld '/' 'tsk' '/']); end
% if ~exist([cfg.clr_fld '/' 'rmv_chn' '/'],'dir'); mkdir([cfg.clr_fld '/' 'rmv_chn' '/']); end
% if ~exist([cfg.clr_fld '/' 'cmn_nse' '/'],'dir'); mkdir([cfg.clr_fld '/' 'cmn_nse' '/']); end
% if ~exist([cfg.clr_fld '/' 'bse_frq' '/'],'dir'); mkdir([cfg.clr_fld '/' 'bse_frq' '/' ]); end
% if ~exist([cfg.clr_fld '/' 'electrode_location' '/'],'dir'); mkdir([cfg.clr_fld '/' 'electrode_location' '/']); end
% if ~exist([cfg.clr_fld '/' 'pial' '/'],'dir'); mkdir([cfg.clr_fld '/' 'pial' '/']); end
% if ~exist([cfg.clr_fld '/' 'hemi' '/'],'dir'); mkdir([cfg.clr_fld '/' 'hemi' '/']); end
% if ~exist([cfg.clr_fld '/' 'trialfun' '/'],'dir'); mkdir([cfg.clr_fld '/' 'trialfun' '/']); end

% if ~exist([cfg.clr_fld '/' 'indir' '/' cfg.sbj_nme],'file'); system(['touch ' cfg.clr_fld '/' 'indir' '/' cfg.sbj_nme]); end
% if ~exist([cfg.clr_fld '/' 'cln_dir' '/' cfg.sbj_nme],'file'); system(['touch ' cfg.clr_fld '/' 'cln_dir' '/' cfg.sbj_nme]); end
% if ~exist([cfg.clr_fld '/' 'end_dir' '/' cfg.sbj_nme],'file'); system(['touch ' cfg.clr_fld '/' 'end_dir' '/' cfg.sbj_nme]); end
% if ~exist([cfg.clr_fld '/' 'tsk' '/' cfg.sbj_nme],'file'); system(['touch ' cfg.clr_fld '/' 'tsk' '/' cfg.sbj_nme]); end
% if ~exist([cfg.clr_fld '/' 'rmv_chn' '/' cfg.sbj_nme],'file'); system(['touch ' cfg.clr_fld '/' 'rmv_chn' '/' cfg.sbj_nme]); end
% if ~exist([cfg.clr_fld '/' 'cmn_nse' '/' cfg.sbj_nme],'file'); system(['touch ' cfg.clr_fld '/' 'cmn_nse' '/' cfg.sbj_nme]); end
% if ~exist([cfg.clr_fld '/' 'bse_frq' '/' cfg.sbj_nme],'file'); system(['touch ' cfg.clr_fld '/' 'bse_frq' '/' cfg.sbj_nme]); end
% if ~exist([cfg.clr_fld '/' 'electrode_location' '/' cfg.sbj_nme],'file'); system(['touch ' cfg.clr_fld '/' 'electrode_location' '/' cfg.sbj_nme]); end
% if ~exist([cfg.clr_fld '/' 'pial' '/' cfg.sbj_nme],'file'); system(['touch ' cfg.clr_fld '/' 'pial' '/' cfg.sbj_nme]); end
% if ~exist([cfg.clr_fld '/' 'hemi' '/' cfg.sbj_nme],'file'); system(['touch ' cfg.clr_fld '/' 'hemi' '/' cfg.sbj_nme]); end
% if ~exist([cfg.clr_fld '/' 'trialfun' '/' cfg.sbj_nme],'file'); system(['touch ' cfg.clr_fld '/' 'trialfun' '/' cfg.sbj_nme]); end

% system(['gedit ' cfg.clr_fld '/' 'indir' '/' cfg.sbj_nme ' &']);
% system(['gedit ' cfg.clr_fld '/' 'cln_dir' '/' cfg.sbj_nme ' &']);
% system(['gedit ' cfg.clr_fld '/' 'end_dir' '/' cfg.sbj_nme ' &']);
% system(['gedit ' cfg.clr_fld '/' 'tsk' '/' cfg.sbj_nme ' &']);
% system(['gedit ' cfg.clr_fld '/' 'rmv_chn' '/' cfg.sbj_nme ' &']);
% system(['gedit ' cfg.clr_fld '/' 'cmn_nse' '/' cfg.sbj_nme ' &']);
% system(['gedit ' cfg.clr_fld '/' 'bse_frq' '/' cfg.sbj_nme ' &']);
% system(['gedit ' cfg.clr_fld '/' 'electrode_location' '/' cfg.sbj_nme ' &']);
% system(['gedit ' cfg.clr_fld '/' 'pial' '/' cfg.sbj_nme ' &']);
% system(['gedit ' cfg.clr_fld '/' 'hemi' '/' cfg.sbj_nme ' &']);