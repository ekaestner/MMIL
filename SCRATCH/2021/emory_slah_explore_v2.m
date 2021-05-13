clear; clc;

dta_dir = '/space/mcdonald-syn01/1/data/MCD_MRI/EmorySLAH';

srf_dir = '/space/mcdonald-syn01/1/data/MCD_MRI/raw';

ejk_chk_dir('/home/ekaestne/Downloads/Emory_Slah_Explore/');

%%
sbj_dir = dir(dta_dir);
sbj_dir = {sbj_dir(:).name};
sbj_dir = sbj_dir(3:end);
sbj_dir( string_find(sbj_dir,'tar.gz')) = [];
sbj_dir( string_find(sbj_dir,'fix.sh')) = [];
sbj_dir( string_find(sbj_dir,'fix_attempts.txt')) = [];
sbj_dir( string_find(sbj_dir,'copy_fix')) = [];

srf_dir_hld = dir( srf_dir );
srf_dte_hld = {srf_dir_hld(:).date};
srf_dir_hld = {srf_dir_hld(:).name};
srf_dir_hld( string_find(srf_dir_hld,'copy_fix')) = [];

%%
sbj_tot_hld = cell(0);
scn_tot_hld = [];
for iS = 1:numel(sbj_dir)
    
    neu_dir = dir( [ dta_dir '/' sbj_dir{iS} ] );
    
    neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})) = {neu_dir(3:end).name};
    
    und_ind = strfind( neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})), '_');
    for iSC = 1:numel(neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})))
        neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})){iSC} = neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})){iSC}(1:und_ind{iSC}(end)-1);
        if strcmpi(neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})){iSC}(end-1:end),'_2'); neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})){iSC} = neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})){iSC}(1:end-2); end
        %         if strcmpi(neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})){iSC}(end-1:end),'_1'); neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})){iSC} = neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})){iSC}(1:end-2); end
        neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})){iSC} = strrep(neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})){iSC},'repeat','Repeat');
        neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})){iSC} = strrep(neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})){iSC},'t1','T1');
    end
    
    scn_tot_hld = [ scn_tot_hld numel(neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'}))) ];
    sbj_tot_hld = { sbj_tot_hld{:} neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})){:} };
    
end

sbj_tot_hld = sort(sbj_tot_hld);
tbl_hld = tabulate(sbj_tot_hld);
sbj_tot_hld = sort(unique(sbj_tot_hld));

cell2csv( '/home/ekaestne/Downloads/Emory_Slah_Explore/emory_slah_table_v2.csv', tbl_hld(:,1:2))

%%
tot_scn_num = max(scn_tot_hld);
tot_sbj_num = max(cell2mat(tbl_hld(:,2)));

scn_out_hld = cell( size(tbl_hld,1), tot_sbj_num+1);
scn_out_hld(:,1) = tbl_hld(:,1);
sbj_out_hld = cell( numel(sbj_dir), tot_scn_num+1);
sbj_out_hld(:,1) = sbj_dir';

for iS = 1:numel(sbj_dir)
    for iSC = 1:numel(neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})))
        
        scn_num = dir( [ dta_dir '/' sbj_dir{iS} '/' neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})){iSC} ] );
        scn_num = num2str(numel( string_find({scn_num(:).name},{'\.dcm$'}) ));
        
        % Subjects
        sbj_out_hld{iS,iSC+1} =[ neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})){iSC} '; ' scn_num] ;
        
        % Scans
        scn_row = find(strcmpi( scn_out_hld(:,1), neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})){iSC}));
        emp_col = find(cellfun(@isempty,scn_out_hld(scn_row,:)),1);
        scn_out_hld{ scn_row, emp_col} = [ sbj_dir{iS} ', ' scn_num];
        
    end
end

cell2csv( '/home/ekaestne/Downloads/Emory_Slah_Explore/emory_slah_scans_v2.csv',    scn_out_hld)
cell2csv( '/home/ekaestne/Downloads/Emory_Slah_Explore/emory_slah_subjects_v2.csv', sbj_out_hld)

%%
scn_hld = mmil_readtext('/home/ekaestne/Downloads/Emory_Slah_Explore/emory_slah_scans_v2.csv');
scn_nme = scn_hld(:,1);

sbj_hld = mmil_readtext('/home/ekaestne/Downloads/Emory_Slah_Explore/emory_slah_subjects_v2.csv');
sbj_col = size(sbj_hld,2);
sbj_col = 2:2:sbj_col;
sbj_nme = sbj_hld(:,1);

%%
out_hld = cell( numel(sbj_nme)+1, numel(scn_nme)+1 );
out_hld(2:end,1) = sbj_nme;
out_hld(1,2:end) = scn_nme;

dte_hld = cell( numel(sbj_nme), 2 );

for iS = 1:numel(sbj_nme)
    
    clear srs_nme_hld csv_nme_hld
    
    sbj_row = find(strcmpi( out_hld(:,1), sbj_nme{iS} ));
    
    srf_ind = string_find( srf_dir_hld, sbj_nme(iS) );
    
    dte_hld{iS,1} = sbj_nme{iS};
    
    if ~isempty(srf_ind)
        
        dte_hld{iS,2} = srf_dte_hld{srf_ind};
        
        srs_inf_mat     = load([ srf_dir '/' srf_dir_hld{srf_ind} '/' 'OrigSeriesInfo.mat']);
        srs_inf_scd_mat = load([ srf_dir '/' srf_dir_hld{srf_ind} '/' 'SeriesInfo.mat']);
        
        for iSI = 1:numel(srs_inf_mat.SeriesInfo)
            srs_nme = srs_inf_mat.SeriesInfo(iSI).FileNames{1};
            sls_ind = strfind(srs_nme,'/');
            srs_nme_hld{iSI} = srs_nme(sls_ind(end-1)+1:sls_ind(end)-1);
            
            und_ind = strfind( srs_nme_hld{iSI}, '_');
            srs_nme_hld{iSI} = srs_nme_hld{iSI}(1:und_ind(end)-1);
            if strcmpi(srs_nme_hld{iSI}(end-1:end),'_2'); srs_nme_hld{iSI} = srs_nme_hld{iSI}(1:end-2); end
            %             if strcmpi(srs_nme_hld{iSI}(end-1:end),'_1'); srs_nme_hld{iSI} = srs_nme_hld{iSI}(1:end-2); end
            srs_nme_hld{iSI} = strrep(srs_nme_hld{iSI},'repeat','Repeat');
            srs_nme_hld{iSI} = strrep(srs_nme_hld{iSI},'t1','T1');
            
            
            srs_nme = srs_inf_scd_mat.SeriesInfo(iSI).FileNames{1};
            sls_ind = strfind(srs_nme,'/');
            csv_nme_hld{iSI} = srs_nme(sls_ind(end-1)+1:sls_ind(end)-1);
        end
        
        if exist([ srf_dir '/' srf_dir_hld{srf_ind} '/' 'SeriesInfo.csv'])
            
            srs_inf_csv = mmil_readtext([ srf_dir '/' srf_dir_hld{srf_ind} '/' 'SeriesInfo.csv']);
            
            for iSC = 2:size( sbj_hld, 2)
                
                scn_typ = sbj_hld{ iS, iSC }(1:strfind(sbj_hld{ iS, iSC },';')-1 );
                
                if ~isempty(scn_typ)
                    scn_col = find(strcmpi( out_hld(1,:), scn_typ));
                    
                    srs_ind = find(strcmpi( srs_nme_hld, scn_typ));
                    if isempty(srs_ind)
                        out_hld{ sbj_row, scn_col } = 'MISSING';
                    else
                        csv_ind = find(strcmpi( srs_inf_csv(:,2),csv_nme_hld{srs_ind(1)}));
                        out_hld{ sbj_row, scn_col } = srs_inf_csv{ csv_ind, 4};
                    end
                end
                
            end
        else
            for iSC = 2:size( sbj_hld, 2)
                scn_typ = sbj_hld{ iS, iSC }(1:strfind(sbj_hld{ iS, iSC },';')-1 );
                if ~isempty(scn_typ)
                    scn_col = find(strcmpi( out_hld(1,:), scn_typ));
                    srs_ind = find(strcmpi( srs_nme_hld, scn_typ));
                    if isempty(srs_ind)
                        out_hld{ sbj_row, scn_col } = 'MISSING';
                    else
                        out_hld{ sbj_row, scn_col } = 'CSV_ERROR';
                    end
                end
            end
        end
    else
        for iSC = 2:size( sbj_hld, 2)
            scn_typ = sbj_hld{ iS, iSC }(1:strfind(sbj_hld{ iS, iSC },';')-1 );
            if ~isempty(scn_typ)
                scn_col = find(strcmpi( out_hld(1,:), scn_typ));
                out_hld{ sbj_row, scn_col } = 'NAME_ERROR';
            end
        end
    end
end

cell2csv( '/home/ekaestne/Downloads/Emory_Slah_Explore/ucsd_preprocessing.csv', out_hld )

tot_out = cell( size(out_hld,1)-1, 3);
for iT = 1:size(out_hld,1)-1
    ind_inc = find(~cellfun(@isempty,out_hld(iT+1,:)));
    nme_inc = strcat(out_hld(iT+1,ind_inc(2:end)),';');
    tot_out{ iT, 2 } = out_hld{iT+1,ind_inc(1)};
    tot_out{ iT, 3 } = cat(2,nme_inc{:});
    
    if isempty( dte_hld{iT,2} )
        tot_out{ iT, 1 } = 'Name error';
    elseif ~isempty( strfind( dte_hld{iT,2}, '21-Jan'))
        tot_out{ iT, 1 } = 'Previously Processed';
    elseif ~isempty( strfind( tot_out{ iT, 3 }, 'UNKNOWN')) || ~isempty( strfind( dte_hld{iT,2}, 'JUNK')) || ~isempty( strfind( dte_hld{iT,2}, 'MISSING'))
        tot_out{ iT, 1 } = 'Unknown Scans';
    elseif ~isempty( strfind( tot_out{ iT, 3 }, 'CSV_ERROR'))
        tot_out{ iT, 1 } = 'CSV Error';
    else
        tot_out{ iT, 1 } = 'Correct?';
    end
end

cell2csv( '/home/ekaestne/Downloads/Emory_Slah_Explore/ucsd_preprocessing_summary.csv', [ dte_hld tot_out ] );

%%
sum_hld = mmil_readtext('/home/ekaestne/Downloads/Emory_Slah_Explore/ucsd_preprocessing_summary.csv');
sbj_hld = mmil_readtext('/home/ekaestne/Downloads/Emory_Slah_Explore/emory_slah_subjects_v2.csv');

sbj_out = cell( size( sbj_hld, 1), 3);
for iS = 1:size( sbj_hld, 1)
    
    sbj_out{iS,1} = sbj_hld{iS,1};
    
    for iSC = 2:size(sbj_hld,2)
        if ~isempty(sbj_hld{iS,iSC})
            scn_add{iSC-1} = sbj_hld{iS,iSC};
            und_ind = strfind( scn_add{iSC-1}, '_');
            scn_add{iSC-1} = scn_add{iSC-1}(1:und_ind(end)-1);
        end
    end
    scn_add = sort(scn_add);
    
    sbj_out{iS,2} = sum_hld{strcmpi(sum_hld(:,1),sbj_out{iS,1}),3};
    
    scn_add = strcat(scn_add,';--');
    scn_add = [scn_add{:}];
    sbj_out{iS,3} = scn_add(1:end-4);
    
    clear scn_add
    
end

cell2csv( '/home/ekaestne/Downloads/Emory_Slah_Explore/emory_slah_subjects_v3.csv', sbj_out)


ttt = tabulate(sbj_out(:,3));
cell2csv( '/home/ekaestne/Downloads/Emory_Slah_Explore/emory_slah_subjects_v3_table.csv', ttt)

%%
ttt = mmil_readtext('/home/ekaestne/Downloads/Emory_Slah_Explore/emory_slah_subjects_v3.csv');

bdd_sbj_ind     = find(strcmpi(ttt(:,2),'CSV Error'));
gdd_sbj_ind     = find(strcmpi(ttt(:,2),'Correct?'));
gdd_sbj_alt_ind = find(strcmpi(ttt(:,2),'Unknown Scans'));

gdd_sbj_hld     = strsplit( ttt{gdd_sbj_ind(7),3}, ';--');
gdd_sbj_alt_hld = strsplit( ttt{gdd_sbj_alt_ind(7),3}, ';--');
err_sbj_hld     = strsplit( ttt{bdd_sbj_ind(7),3}, ';--');
setxor( gdd_sbj_alt_hld, err_sbj_hld)

%% Header error investigation
sbj_hld = mmil_readtext('/home/ekaestne/Downloads/Emory_Slah_Explore/emory_slah_subjects_v3.csv');

sbj_dir = dir(srf_dir);
sbj_dir = {sbj_dir(:).name};

for iS = 1:size(sbj_hld,1)
    
    sbj_dir_ind = string_find( sbj_dir, sbj_hld{iS,1});
    
    scn_dir = dir([ srf_dir '/' sbj_dir{sbj_dir_ind}]);
    scn_dir = {scn_dir(:).name};
    scn_dir = scn_dir( string_find( scn_dir, {'st001'}) );
    
    if ~isempty(scn_dir)
        dcm_fle = [ srf_dir '/' sbj_dir{sbj_dir_ind} '/' scn_dir{1} '/' 'im0001.dcm'];
        
        cmd = [ 'dcmdump ' dcm_fle ' | ' 'grep ''DateOfLastCalibration''' ];
        
        [ ~, out_put{iS}] = unix(cmd);
        out_put{iS} = strrep( out_put{iS}, ',', ';');
    end
end

cell2csv( '/home/ekaestne/Downloads/Emory_Slah_Explore/header_investigation.csv', [ sbj_hld(:,1:2) out_put']);
clear out_put

%% Try to fix the header problem
usd_pre_prc = mmil_readtext('/home/ekaestne/Downloads/Emory_Slah_Explore/ucsd_preprocessing_summary.csv');

dta_dir = '/space/mcdonald-syn01/1/data/MCD_MRI/EmorySLAH/'; %FAMOUSFACESEPK183_1102_BITC20160825_copy_fix_attempt';

for iS = 1:size(usd_pre_prc,1)
    
    if strcmpi( usd_pre_prc{iS,3}, 'CSV Error')
        
        dta_dir_sbj = [ dta_dir '/' usd_pre_prc{iS,1} ];
        
        dta_img_dir = dir([ dta_dir_sbj '/']);
        dta_img_dir = {dta_img_dir(:).name};
        dta_img_dir = dta_img_dir(3:end);
        
        dta_dir_sbj = [ dta_dir_sbj '_' 'copy_fix' ];
        ejk_chk_dir(dta_dir_sbj);
        
        for iF = 1:numel(dta_img_dir)
            
            if isempty( strfind( dta_img_dir{iF}, '_ColFA') ) && ...
                    isempty( strfind( dta_img_dir{iF}, '_ADC') ) && ...
                    isempty( strfind( dta_img_dir{iF}, '_FA') ) && ...
                    isempty( strfind( dta_img_dir{iF}, '_TENSOR') ) && ...
                    isempty( strfind( dta_img_dir{iF}, '_TRACEW') )
                
                unix( [ 'cp -r ' dta_dir '/' usd_pre_prc{iS,1} '/' dta_img_dir{iF} ' ' dta_dir_sbj '/' dta_img_dir{iF} '/' ] );
                
                dta_img_fle = dir( [ dta_dir_sbj '/' dta_img_dir{iF} ] );
                dta_img_fle = {dta_img_fle(:).name};
                dta_img_fle = dta_img_fle( ~cellfun(@isempty,regexp(dta_img_fle,'.dcm$')));
                
                cd([ dta_dir_sbj '/' dta_img_dir{iF} ]);
                
                %     if 0
                
                tic
                for iD = 1:numel(dta_img_fle)
                    
                    cmd_hld = 'dcmodify -m';
                    
                    cmd_run= [ cmd_hld ' "(0018,1200)=20090304" "' dta_img_fle{iD} '" -nb -imt -q'];
                    unix(cmd_run);
                    
                    cmd_run = [ cmd_hld ' "(0018,1201)=123723.000000" "' dta_img_fle{iD} '" -nb -imt -q'];
                    [ttt, tt2] = unix(cmd_run);
                    
                end
                
                tme_hld = toc;
                fprintf('Subject %i of %i, scan %i of %i - %s : %i\n', iS, size(usd_pre_prc,1), iF, numel(dta_img_dir), usd_pre_prc{iS,1},tme_hld)
                
                %     else
                %
                %         fprintf('%s : skipping\n',[ dta_dir '/' dta_img_dir{iF} ])
                %
                %     end
                
            end
        end
    end
end

%% Quick open for Don
srs_inf_dir = '/space/mcdonald-syn01/1/data/MCD_MRI/raw/';
srs_inf_fld = dir(srs_inf_dir);
srs_inf_fld = {srs_inf_fld(:).name};
dbs_out_cmr_dir = '/space/mcdonald-syn01/1/data/MCD_MRI/batchdirs/MCD_MRI_preproc_sMRI_orig_new_20210429_160415/pbsout/';
dbs_out_mtc_dir = '/space/mcdonald-syn01/1/data/MCD_MRI/batchdirs/MCD_MRI_preproc_sMRI_orig_new_20210422_133912/pbsout/';



sbj_nme = 'FAMOUSFACESEPK238_1205_BITC20160818';

srs_inf_ind     = string_find( srs_inf_fld, sbj_nme );
unix(['oocalc ' srs_inf_dir '/' srs_inf_fld{srs_inf_ind} '/' 'SeriesInfo.csv']);

[ ~, dbs_out_cmr_ind] = unix(['grep -R ''' sbj_nme '''' ' ' dbs_out_cmr_dir ]);
[ ~, dbs_out_mtc_ind] = unix(['grep -R ''' sbj_nme '''' ' ' dbs_out_mtc_dir ]);
if isempty(dbs_out_cmr_ind)
    ttt = strfind( dbs_out_mtc_ind, 'job_' );
    sss = strfind( dbs_out_mtc_ind, '.out' );
    unix(['kwrite ' dbs_out_mtc_dir '/' dbs_out_mtc_ind(ttt(1):sss(1)+3) ]);
else
    ttt = strfind( dbs_out_cmr_ind, 'job_' );
    sss = strfind( dbs_out_cmr_ind, '.out' );
    unix(['kwrite ' dbs_out_cmr_dir '/' dbs_out_cmr_ind(ttt(1):sss(1)+3) ]);
end

%% Final spreadsheet (??)
sbj_dir = dir(dta_dir);
sbj_dir = {sbj_dir(:).name};

scn_cls = mmil_readtext('/home/ekaestne/Downloads/Emory_Slah_Explore/ucsd_preprocessing.csv');
scn_cls_col = scn_cls(1,:); scn_cls_col{1} = '';
scn_cls     = scn_cls(2:end,:);

sum_scn_col = { 'bitc_diff__Matt_LT_LB_60dn'       'bitc_diff__Matt_LT_LB_60up' ...
    'cmrr_mb_diff1.75mmMB3AP'          'cmrr_mb_diff1.75mmMB3PA' ...
    'cmrr_mbep2d_diff1.75mmAP3Bvalues' 'cmrr_mbep2d_diff1.75mmPA3Bvalues' };
sum_scn_col_ind = find(ismember( scn_cls_col, sum_scn_col ));

scn_seq_nme = { 'bitc' '' 'cmrr_mb' '' 'cmrr_mbep2d' '' };

col_nme = { 'SubjID' 'Sequence' 'For_Outcome' 'Rev_Outcome' ...
    'For_Name' 'For_Num' 'For_B0' 'For_B1000' 'For_B2000' ...
    'Rev_Name' 'Rev_Num' 'Rev_B0' 'Rev_B1000' 'Rev_B2000'};
out_tbl = cell( size(scn_cls,1), numel(col_nme));
for iS = 48:size(scn_cls,1)
    
    sbj_dir_use = sbj_dir( string_find( sbj_dir, scn_cls{iS,1}) );
    sbj_dir_use( string_find(sbj_dir_use,'copy_fix') ) = [];
    
    dti_scn_col = find( ~cellfun(@isempty,scn_cls( iS, sum_scn_col_ind)) );
    if numel(dti_scn_col)==2
        
        scn_dir = dir([ dta_dir '/' sbj_dir_use{1}]);
        scn_dir = {scn_dir(:).name};
        
        % FIRST %%%%%%
        out_tbl{iS,1} = scn_cls{iS,1};
        
        out_tbl{iS,2} = scn_seq_nme{ dti_scn_col(1) };
        
        out_tbl{iS,3} = scn_cls{ iS, sum_scn_col_ind(dti_scn_col(1)) };
        out_tbl{iS,4} = scn_cls{ iS, sum_scn_col_ind(dti_scn_col(2)) };
        %%%%%%%
        
        % FORWARD %%%%%%
        scn_dir_for = find(~cellfun(@isempty,regexp( scn_dir, [ '^' sum_scn_col{dti_scn_col(1)} '_\d'] )));
        
        tmp_dir_hld = dir( [ dta_dir '/' sbj_dir_use{1} '/' scn_dir{scn_dir_for} '/'] );
        tmp_dir_hld = {tmp_dir_hld(:).name};
        tmp_dir_hld = tmp_dir_hld(string_find(tmp_dir_hld,'.dcm$'));
        
        dcm_inf = dicominfo([ dta_dir '/' sbj_dir_use{1} '/' scn_dir{scn_dir_for} '/' tmp_dir_hld{1}]);
        out_tbl{iS,5} = dcm_inf.ProtocolName;
        
        out_tbl{iS,6} = numel(tmp_dir_hld);
        
        [ ~, dcm_hed] = unix(['dcmdump ' dta_dir '/' sbj_dir_use{1} '/' scn_dir{scn_dir_for} '/' '*.dcm | grep -i ''0019,100c'' | sort']);
        dcm_hed_beg = strfind(dcm_hed,'IS [');
        dcm_hed_end = regexp(dcm_hed,'\d]');
        bvl_hld = nan(numel(dcm_hed_beg),1);
        for iN = 1:numel(dcm_hed_beg)
            bvl_hld(iN) = str2double(dcm_hed(dcm_hed_beg(iN)+4:dcm_hed_end(iN)));
        end
        out_tbl{iS,7} = sum(bvl_hld<50);
        out_tbl{iS,8} = sum(bvl_hld>950 & bvl_hld<1050);
        out_tbl{iS,9} = sum(bvl_hld>1950 & bvl_hld<2050);
        %%%%%%%
        
        % REVERSE %%%%%%
        scn_dir_rev = find(~cellfun(@isempty,regexp( scn_dir, [ '^' sum_scn_col{dti_scn_col(2)} '_\d'] )));
        
        tmp_dir_hld = dir( [ dta_dir '/' sbj_dir_use{1} '/' scn_dir{scn_dir_rev} '/'] );
        tmp_dir_hld = {tmp_dir_hld(:).name};
        tmp_dir_hld = tmp_dir_hld(string_find(tmp_dir_hld,'.dcm$'));
        
        dcm_inf = dicominfo([ dta_dir '/' sbj_dir_use{1} '/' scn_dir{scn_dir_rev} '/' tmp_dir_hld{1}]);
        out_tbl{iS,10} = dcm_inf.ProtocolName;
        
        out_tbl{iS,11} = numel(tmp_dir_hld);
        
        [ ~, dcm_hed] = unix(['dcmdump ' dta_dir '/' sbj_dir_use{1} '/' scn_dir{scn_dir_rev} '/' '*.dcm | grep -i ''0019,100c'' | sort']);
        dcm_hed_beg = strfind(dcm_hed,'IS [');
        dcm_hed_end = regexp(dcm_hed,'\d]');
        bvl_hld = nan(numel(dcm_hed_beg),1);
        for iN = 1:numel(dcm_hed_beg)
            bvl_hld(iN) = str2double(dcm_hed(dcm_hed_beg(iN)+4:dcm_hed_end(iN)));
        end
        out_tbl{iS,12} = sum(bvl_hld<50);
        out_tbl{iS,13} = sum(bvl_hld>950 & bvl_hld<1050);
        out_tbl{iS,14} = sum(bvl_hld>1950 & bvl_hld<2050);
        %%%%%%%
    elseif numel(dti_scn_col==1)
        dti_scn_col = [dti_scn_col-1 dti_scn_col];
        
        out_tbl{iS,1} = scn_cls{iS,1};
        
        out_tbl{iS,2} = scn_seq_nme{ dti_scn_col(1) };
        
        out_tbl{iS,3} = 'MISSING';
        out_tbl{iS,4} = scn_cls{ iS, sum_scn_col_ind(dti_scn_col(2)) };
    else
        out_tbl{iS,1} = scn_cls{iS,1};
        out_tbl{iS,2} = 'MISSING';
        out_tbl{iS,3} = 'MISSING';
        out_tbl{iS,4} = 'MISSING';
    end
end

cell2csv('/home/ekaestne/Downloads/Emory_Slah_Explore/Emory_Scan_Classifications.csv', [ col_nme ; out_tbl ]);

%% Classify add









