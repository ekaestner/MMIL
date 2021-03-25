function iSASZ_stat_table(fcfg)

% fcfg.typ
%
% fcfg.sig_in__pth
% fcfg.tab_out_pth
% fcfg.tab_fle_nme
%
% fcfg.tab_col
% fcfg.ovr_lap_col
% fcfg.ovr_lap_lbl
%
% fcfg.plt_out_pth
% fcfg.plt_fle_nme
% fcfg.cmb_plt_typ
% fcfg.plt_grp
% fcfg.plt_grp_nme
%
% fcfg.loc_in__pth

if ~isfield(fcfg,'str_sbj'); fcfg.str_sbj = 1; end

sbj_nme = dir(fcfg.sig_in__pth); sbj_nme = {sbj_nme([sbj_nme.isdir]).name}; sbj_nme = sbj_nme(3:end);

for iRM = 1:numel(fcfg.plt_grp_nme)
    %if isdir([fcfg.plt_out_pth '/' fcfg.plt_grp_nme{iRM}]); rmdir([fcfg.plt_out_pth '/' fcfg.plt_grp_nme{iRM}],'s'); end
    mkdir([fcfg.plt_out_pth '/' fcfg.plt_grp_nme{iRM}])
end

col_nme = [fcfg.tab_col fcfg.ovr_lap_lbl];

%% Loop ECOG through tables
for iD = 1:numel(fcfg.typ)
    
    if fcfg.ecg == 1;
        
        cur_tab = ['subject'  fcfg.tab_col fcfg.ovr_lap_lbl 'Total Channels'];
        
        for iCL = 1:numel(col_nme)
            ele_loc_hld.(col_nme{iCL}) = [];
        end
        tot_ele_loc = [];        
        
        for iT = fcfg.str_sbj:numel(sbj_nme);
            
            for iTE = 1:numel(fcfg.plt_grp_nme)
                ele_loc_tab.(fcfg.plt_grp_nme{iTE}).(sbj_nme{iT}) = [];
            end
            
            tab_fle = dir([fcfg.sig_in__pth '/' sbj_nme{iT} '/' ]); tab_fle = {tab_fle(:).name}; tab_fle = tab_fle(3:end); tab_fle = tab_fle(~cellfun(@isempty,strfind(tab_fle,'ecog')));
            tab_fle = mmil_readtext([fcfg.sig_in__pth '/' sbj_nme{iT} '/' tab_fle{find(~cellfun('isempty',strfind(tab_fle,fcfg.typ{iD})))}]);
            
            cur_tab{iT+1,1} = sbj_nme{iT};
            
            % Get counts
            for iC = 1:numel(fcfg.tab_col)
                
                if sum(strcmpi(fcfg.tab_col{iC},tab_fle(1,:))) > 0
                    cur_tab{iT+1,iC+1} = sum(cell2mat(tab_fle(2:end,find(strcmpi(fcfg.tab_col{iC},tab_fle(1,:))))));
                    if isfield(fcfg,'ovr_lap_col')
                        ovr_lap_hld(:,iC) = cell2mat(tab_fle(2:end,find(strcmpi(fcfg.tab_col{iC},tab_fle(1,:)))));
                    end
                else
                    cur_tab{iT+1,iC+1} = '-';
                    if isfield(fcfg,'ovr_lap_col')
                        ovr_lap_hld(:,iC) = repmat(-1,size(tab_fle,1)-1,1);
                    end
                end
                
            end
            
            % Get overlaps
            for iO = 1:numel(fcfg.ovr_lap_lbl)
                
                if ~any(any(ovr_lap_hld(:,fcfg.ovr_lap_col{iO}) == -1))
                    cur_tab{iT+1,iC+iO+1} = sum(sum(ovr_lap_hld(:,fcfg.ovr_lap_col{iO}),2)==numel(fcfg.ovr_lap_col{iO}));
                else
                    cur_tab{iT+1,iC+iO+1} = '-';
                end
                
            end
            
            for iO = 1:numel(fcfg.ovr_lap_lbl)
                for iCH = 1:size(ovr_lap_hld,1)
                    ovr_lap_hld(iCH,iC+iO) = ovr_lap_hld(iCH,fcfg.ovr_lap_col{iO}(1)) & ovr_lap_hld(iCH,fcfg.ovr_lap_col{iO}(2));
                end                
            end
            
            cur_tab{iT+1,end} = size(tab_fle,1);
            
            % Record Keeping for plotting
            plt_hld = cell(1,numel(fcfg.plt_grp_nme));
            for iR = 1:numel(fcfg.plt_grp)
                
                plt_num_hld = [];
                for iRC = 1:numel(fcfg.plt_grp{iR})
                    plt_num_hld = [plt_num_hld ; tab_fle(find(ovr_lap_hld(:,fcfg.plt_grp{iR}(iRC))==1)+1,2)];
                end
                plt_num_hld = unique(plt_num_hld);
                
                plt_hld{iR} = [plt_hld{iR} ; plt_num_hld];
                
            end
                        
            % Make electrode location table holder
            loc_fle = dir([fcfg.loc_in__pth '/'  'electrode_anatomical_location' '/' ]); loc_fle = {loc_fle(:).name}; loc_fle = loc_fle(~cellfun(@isempty,strfind(loc_fle,'ecog'))); loc_fle = loc_fle(~cellfun(@isempty,strfind(loc_fle,sbj_nme{iT})));

            if ~isempty(loc_fle) 
                
                if numel(loc_fle) > 1
                    
                    loc_fle1 = mmil_readtext([fcfg.loc_in__pth '/'  'electrode_anatomical_location' '/' loc_fle{1}],'[ ,]');
                    loc_fle2 = mmil_readtext([fcfg.loc_in__pth '/'  'electrode_anatomical_location' '/' loc_fle{2}],'[ ,]');
                    if size(loc_fle1,2) < size(loc_fle2,2)
                        loc_fle = [loc_fle1 cell(size(loc_fle1,1),size(loc_fle2,2)-size(loc_fle1,2)); loc_fle2];
                    elseif size(loc_fle1,2) > size(loc_fle2,2)
                        loc_fle = [loc_fle1 ; loc_fle2 cell(size(loc_fle2,1),size(loc_fle1,2)-size(loc_fle2,2))];
                    else
                        loc_fle = [loc_fle1 ; loc_fle2];
                    end
                 else
                    loc_fle = mmil_readtext([fcfg.loc_in__pth '/'  'electrode_anatomical_location' '/' loc_fle{1}],'[ ,]');
                end
                
                for iCL = 1:size(loc_fle,1)
                    ele_loc(iCL,1) = loc_fle(iCL,1);
                    if ~strcmpi(loc_fle(iCL,5),'100.00%') && ~isempty(loc_fle{iCL,5})
                        loc_pct = loc_fle(iCL,[5:2:end]); loc_pct(cellfun(@isempty,loc_pct)) = [];
                        loc_pct = cellfun(@str2num,cellfun(@(x,y) x(1:y-1),loc_pct,strfind(loc_pct,'.'),'uni',0));
                        [~,loc_pct] = max(loc_pct);
                        ele_loc(iCL,2) = loc_fle(iCL,4+loc_pct*2);
                    else
                        ele_loc(iCL,2) = loc_fle(iCL,6);
                    end                       
                end
                
                for iPL = 1:numel(col_nme)
                    for iCL = 1:size(ovr_lap_hld,1)
                        ele_loc_ind = find(strcmpi(ele_loc(:,1),tab_fle(iCL+1,2)));
                        if ovr_lap_hld(iCL,iPL); ele_loc_hld.(col_nme{iPL}) = [ele_loc_hld.(col_nme{iPL}) ; ele_loc(ele_loc_ind,2)]; end
                    end
                end
                
                % Record Keeping for location plotting
                for iCN = 1:numel(col_nme)
                    loc_hld.(sbj_nme{iT}).(col_nme{iCN}) = [];
                end
                
                for iLC = 1:numel(col_nme)
                    if ovr_lap_hld(1,iLC) >-1
                        loc_hld.(sbj_nme{iT}).(col_nme{iLC}) = tab_fle(find(ovr_lap_hld(:,iLC)==1)+1,2);
                    else
                        loc_hld.(sbj_nme{iT}) = rmfield(loc_hld.(sbj_nme{iT}),col_nme{iLC});
                    end
                end
                
                tot_ele_loc = [tot_ele_loc ele_loc(:,2)'];
                
            end
            
            % Make electrode plot mover holder
            for iTE = 1:numel(fcfg.plt_grp_nme)
                ele_loc_tab.(fcfg.plt_grp_nme{iTE}).(sbj_nme{iT}) = [ele_loc_tab.(fcfg.plt_grp_nme{iTE}).(sbj_nme{iT}) ; [repmat({sbj_nme{iT}},numel(plt_hld{iTE}),1) plt_hld{iTE}]];
            end            

            clear ovr_lap_hld
            
        end        
        
        % Make Anatomical Table
        ana_nme = unique(unique(tot_ele_loc(~cellfun(@isempty,tot_ele_loc))));
        ana_tbl = {'Location' col_nme{:} 'Total'};
        for iA = 1:numel(ana_nme)
            ana_tbl{iA+1,1} = ana_nme{iA};
            ana_tbl{iA+1,numel(col_nme)+2} = sum(strcmpi(tot_ele_loc,ana_nme{iA}));
            for iE = 1:numel(col_nme)
                ana_tbl{iA+1,iE+1} = sum(strcmpi(ele_loc_hld.(col_nme{iE}),ana_nme{iA}))/ana_tbl{iA+1,numel(col_nme)+2};
            end            
        end
        
        clear tot_ele_loc
        
        % save table
        cell2csv([fcfg.tab_out_pth '/' fcfg.tab_fle_nme '_' fcfg.typ{iD} '_ecog.csv'],cur_tab)
        cell2csv([fcfg.tab_out_pth '/' fcfg.tab_fle_nme '_' fcfg.typ{iD} '_ecog_location.csv'],ana_tbl)
        save([fcfg.tab_out_pth '/' fcfg.tab_fle_nme '_' fcfg.typ{iD} '_ecog_plot.mat'],'ele_loc_tab')
        save([fcfg.tab_out_pth '/' fcfg.tab_fle_nme '_' fcfg.typ{iD} '_ecog_location_plot.mat'],'loc_hld')
        
    end
    
end

%% Loop Depth through tables
for iD = 1:numel(fcfg.typ)
    
    if fcfg.dep == 1;
        
        cur_tab = ['subject'  fcfg.tab_col fcfg.ovr_lap_lbl 'Total Channels'];
        
        for iCL = 1:numel(col_nme)
            ele_loc_hld.(col_nme{iCL}) = [];
        end
        tot_ele_loc = [];        
        
        for iT = fcfg.str_sbj:numel(sbj_nme);
            
            for iTE = 1:numel(fcfg.plt_grp_nme)
                ele_loc_tab.(fcfg.plt_grp_nme{iTE}).(sbj_nme{iT}) = [];
            end
            
            tab_fle = dir([fcfg.sig_in__pth '/' sbj_nme{iT} '/' ]); tab_fle = {tab_fle(:).name}; tab_fle = tab_fle(3:end); tab_fle = tab_fle(~cellfun(@isempty,strfind(tab_fle,'depth')));
            tab_fle = mmil_readtext([fcfg.sig_in__pth '/' sbj_nme{iT} '/' tab_fle{find(~cellfun('isempty',strfind(tab_fle,fcfg.typ{iD})))}]);
            
            cur_tab{iT+1,1} = sbj_nme{iT};
            
            if size(tab_fle,1) > 1
            
            % Get counts
            for iC = 1:numel(fcfg.tab_col)
                
                if sum(strcmpi(fcfg.tab_col{iC},tab_fle(1,:))) > 0
                    cur_tab{iT+1,iC+1} = sum(cell2mat(tab_fle(2:end,find(strcmpi(fcfg.tab_col{iC},tab_fle(1,:))))));
                    if isfield(fcfg,'ovr_lap_col')
                        ovr_lap_hld(:,iC) = cell2mat(tab_fle(2:end,find(strcmpi(fcfg.tab_col{iC},tab_fle(1,:)))));
                    end
                else
                    cur_tab{iT+1,iC+1} = '-';
                    if isfield(fcfg,'ovr_lap_col')
                        ovr_lap_hld(:,iC) = repmat(-1,size(tab_fle,1)-1,1);
                    end
                end
                
            end
            
            % Get overlaps
            for iO = 1:numel(fcfg.ovr_lap_lbl)
                
                if ~any(any(ovr_lap_hld(:,fcfg.ovr_lap_col{iO}) == -1))
                    cur_tab{iT+1,iC+iO+1} = sum(sum(ovr_lap_hld(:,fcfg.ovr_lap_col{iO}),2)==numel(fcfg.ovr_lap_col{iO}));
                else
                    cur_tab{iT+1,iC+iO+1} = '-';
                end
                
            end
            
            for iO = 1:numel(fcfg.ovr_lap_lbl)
                for iCH = 1:size(ovr_lap_hld,1)
                    ovr_lap_hld(iCH,iC+iO) = ovr_lap_hld(iCH,fcfg.ovr_lap_col{iO}(1)) & ovr_lap_hld(iCH,fcfg.ovr_lap_col{iO}(2));
                end                
            end
            
            cur_tab{iT+1,end} = size(tab_fle,1);
            
            % Record Keeping for plotting
            plt_hld = cell(1,numel(fcfg.plt_grp_nme));
            for iR = 1:numel(fcfg.plt_grp)
                
                plt_num_hld = [];
                for iRC = 1:numel(fcfg.plt_grp{iR})
                    plt_num_hld = [plt_num_hld ; tab_fle(find(ovr_lap_hld(:,fcfg.plt_grp{iR}(iRC))==1)+1,2)];
                end
                plt_num_hld = unique(plt_num_hld);
                
                plt_hld{iR} = [plt_hld{iR} ; plt_num_hld];
                
            end
                                    
            % Make electrode location table holder
            loc_fle = dir([fcfg.loc_in__pth '/'  'electrode_anatomical_location' '/' ]); loc_fle = {loc_fle(:).name}; loc_fle = loc_fle(~cellfun(@isempty,strfind(loc_fle,'depth'))); loc_fle = loc_fle(~cellfun(@isempty,strfind(loc_fle,sbj_nme{iT})));

            if ~isempty(loc_fle) 
                
                if numel(loc_fle) > 1
                    try loc_fle1 = mmil_readtext([fcfg.loc_in__pth '/'  'electrode_anatomical_location' '/' loc_fle{1}],'[ ,]'); catch loc_fle1 = []; end
                    try loc_fle2 = mmil_readtext([fcfg.loc_in__pth '/'  'electrode_anatomical_location' '/' loc_fle{2}],'[ ,]'); catch loc_fle2 = []; end
                    if size(loc_fle1,2) < size(loc_fle2,2)
                        loc_fle = [loc_fle1 cell(size(loc_fle1,1),size(loc_fle2,2)-size(loc_fle1,2)); loc_fle2];
                    elseif size(loc_fle1,2) > size(loc_fle2,2)
                        loc_fle = [loc_fle1 ; loc_fle2 cell(size(loc_fle2,1),size(loc_fle1,2)-size(loc_fle2,2))];
                    else
                        loc_fle = [loc_fle1 ; loc_fle2];
                    end
                else
                    loc_fle = mmil_readtext([fcfg.loc_in__pth '/'  'electrode_anatomical_location' '/' loc_fle{1}],'[ ,]');
                end
                
                for iCL = 1:size(loc_fle,1)
                    ele_loc(iCL,1) = loc_fle(iCL,1);
                    if ~strcmpi(loc_fle(iCL,5),'100.00%') && ~isempty(loc_fle{iCL,5})
                        loc_pct = loc_fle(iCL,[5:2:end]); loc_pct(cellfun(@isempty,loc_pct)) = [];
                        loc_pct = cellfun(@str2num,cellfun(@(x,y) x(1:y-1),loc_pct,strfind(loc_pct,'.'),'uni',0));
                        [~,loc_pct] = max(loc_pct);
                        ele_loc(iCL,2) = loc_fle(iCL,4+loc_pct*2);
                    else
                        ele_loc(iCL,2) = loc_fle(iCL,6);
                    end                       
                end
                
                for iPL = 1:numel(col_nme)
                    for iCL = 1:size(ovr_lap_hld,1)
                        ele_loc_ind = find(strcmpi(ele_loc(:,1),tab_fle(iCL+1,2)));
                        if ovr_lap_hld(iCL,iPL); ele_loc_hld.(col_nme{iPL}) = [ele_loc_hld.(col_nme{iPL}) ; ele_loc(ele_loc_ind,2)]; end
                    end
                end
                
                tot_ele_loc = [tot_ele_loc ele_loc(:,2)'];
                
            end
            
            % Make electrode plot mover holder
            for iTE = 1:numel(fcfg.plt_grp_nme)
                ele_loc_tab.(fcfg.plt_grp_nme{iTE}).(sbj_nme{iT}) = [ele_loc_tab.(fcfg.plt_grp_nme{iTE}).(sbj_nme{iT}) ; [repmat({sbj_nme{iT}},numel(plt_hld{iTE}),1) plt_hld{iTE}]];
            end            
            
            end
            
            clear ovr_lap_hld
            
        end        
        
        % Make Anatomical Table
        ana_nme = unique(unique(tot_ele_loc(~cellfun(@isempty,tot_ele_loc))));
        ana_tbl = {'Location' col_nme{:} 'Total'};
        for iA = 1:numel(ana_nme)
            ana_tbl{iA+1,1} = ana_nme{iA};
            ana_tbl{iA+1,numel(col_nme)+2} = sum(strcmpi(tot_ele_loc,ana_nme{iA}));
            for iE = 1:numel(col_nme)
                ana_tbl{iA+1,iE+1} = sum(strcmpi(ele_loc_hld.(col_nme{iE}),ana_nme{iA}))/ana_tbl{iA+1,numel(col_nme)+2};
            end            
        end
        
        clear tot_ele_loc
        
        % save table
        cell2csv([fcfg.tab_out_pth '/' fcfg.tab_fle_nme '_' fcfg.typ{iD} '_depth.csv'],cur_tab)
        cell2csv([fcfg.tab_out_pth '/' fcfg.tab_fle_nme '_' fcfg.typ{iD} '_depth_location.csv'],ana_tbl)
        save([fcfg.tab_out_pth '/' fcfg.tab_fle_nme '_' fcfg.typ{iD} '_depth_plot.mat'],'ele_loc_tab')
        
    end
    
end

end