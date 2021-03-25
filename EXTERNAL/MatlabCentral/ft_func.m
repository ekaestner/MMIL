%% Function to create a nested dynamic fieldnames structure and pass it to ft_func_v2
% fcfg.data_name = 'all' (default), or specific fields (for 2 datasets at the
%   same time, must be in a n-x-dataset cell)
% fcfg.data_new = 'yes' or 'no' (default) tells whether to overwrite old structure (if
%   filtering etc.) or to create a new structure (if calculating HGP, etc.)
% fcfg.multi    = cell with two rows. 1st identifies which fields are cells (to allow unique
%   arguments to be given to specific data structures, 2nd identifies order
%   of multiple fields |example: {'dataname' 'trl';[1 2 1 2] [1 2 3 4]}|
% fcfg.rmfield = remove fields 'yes'/'no'(default)
% fcfg.empty   = whether or not to try and get a return value ('no' or default'yes')
% fcfg.new_name = new name for the data to be used (only for use with
%   fcfg.new = 'yes')
% fcfg.new_suffix = new suffix to be added to data (only for use with
%   fcfg.new = 'yes')
%
%
% Created by Erik Kaestner (ekaestne@ucsd.edu) 4/3/14

function varargout = ft_func(func,fcfg,dat)

cfg = fcfg;

%% FIND OUT WHICH DATASETS TO WORK WITH
if ~isfield(fcfg,'cmb_dat') && ~isfield(fcfg,'save') && (~isfield(fcfg,'load') || ~strcmpi(fcfg.load,'yes'))
    if ~isfield(cfg,'data_name') && isvar('dat') || isfield(cfg,'data_name') && any(strcmpi(cfg.data_name,'all'));
        cfg.data_name = dat.data_name;
        dat_ind = cellfun(@find,cellfun(@(x) strcmpi([dat.data_name],x),cfg.data_name,'uni',0));
    elseif isfield(cfg,'data_name') && isnumeric(cfg.data_name) % If data_name is a numeric indice
        dat_ind = cfg.data_name;
    else % In special cases where data is not sent in (currenlty only ft_definetrial)
        dat_ind = 1:numel(cfg.dataset);
    end
    if ~isfield(cfg,'data_new'); cfg.data_new = []; end
elseif isfield(fcfg,'cmb_dat') %% Combine different datasets for plotting, stats, etc.
    cfg.rtn_trl    = 0;
    trl_off(1) = 0;
    trl_off(2) = numel(dat.(dat.data_name{fcfg.cmb_dat(1)}).trialinfo);
    if isfield(fcfg,'plt_alt_eve'); dat.(dat.data_name{fcfg.cmb_dat(1)}).trialinfo = dat.(dat.data_name{fcfg.cmb_dat(1)}).cfg.alt_eve.(fcfg.plt_alt_eve); end
    for iCM  = fcfg.cmb_dat(2:end)
        if isfield(fcfg,'plt_alt_eve'); dat.(dat.data_name{fcfg.cmb_dat(iCM)}).trialinfo = dat.(dat.data_name{fcfg.cmb_dat(iCM)}).cfg.alt_eve.(fcfg.plt_alt_eve); end
        ccfg = [];
        ccfg.plt_hck       = [];
        ccfg.data_name     = iCM;
        ccfg.return_events = 0;
        ccfg.old_events    = num2cell(unique(dat.(dat.data_name{iCM}).trialinfo)');
        ccfg.new_events    = num2cell(unique(dat.(dat.data_name{iCM}).trialinfo)'+1000*(iCM-1));
        dat = ft_func(@ft_redefine_events,ccfg,dat);
        trl_off(iCM+1) = numel(dat.(dat.data_name{fcfg.cmb_dat(iCM)}).trialinfo);
    end  
    
    ccfg           = [];
    ccfg.data_name = [ones(1,numel(fcfg.cmb_dat)-1)*fcfg.cmb_dat(1) ; ones(1,numel(fcfg.cmb_dat)-1).*fcfg.cmb_dat(2:end)];
    ccfg.data_new  = 'yes';
    ccfg.methapp   = 'channels';
    if isfield(dat.(dat.data_name{fcfg.cmb_dat(1)}),'trial'); dat = ft_func(@ft_appenddata,ccfg,dat); elseif isfield(dat.(dat.data_name{fcfg.cmb_dat(1)}),'fourierspctrm'); ccfg.parameter = 'fourierspctrm'; dat = ft_func(@ft_appendfreq,ccfg,dat); end
    if isfield(fcfg,'plt_alt_eve'); dat.(dat.data_name{fcfg.cmb_dat(1)}).cfg.alt_eve.(fcfg.plt_alt_eve) = dat.(dat.data_name{fcfg.cmb_dat(1)}).trialinfo; end  
    
    dat_ind = fcfg.cmb_dat(1);
    
    for iCN = 1:size(cfg.cmb_idn,2);
        for iRW = 1:size(cfg.cmb_idn,1);
            cfg.events{iRW,iCN} = cfg.events{iRW,iCN}+1000*(cfg.cmb_idn{iRW,iCN}-1);
            cfg.trls{iRW,iCN}   = cfg.trls{iRW,iCN} + sum(trl_off(1:iCN));
        end
    end
    
    if ~isfield(cfg,'data_new'); cfg.data_new = []; end
    cfg.rtn_trl = 1;
end

%% WORK WITH DATA
if isfield(cfg,'rmfield') && strcmpi(cfg.rmfield,'yes') % Throw away cells no longer needed
    
    varargout{1} = dat;
    
    if isfield(cfg,'sve_fld')
        for iRM = 1:size(dat_ind,2)
            for iSV = 1:numel(cfg.sve_fld)
                varargout{1}.sve_fld.([varargout{1}.data_name{dat_ind(iRM)} '_' cfg.sve_fld{iSV}]) = varargout{1}.(varargout{1}.data_name{dat_ind(iRM)}).cfg.(cfg.sve_fld{iSV});
            end
       end
    end
    
    varargout{1} = rmfield(varargout{1},varargout{1}.data_name(dat_ind));
    varargout{1}.data_name(dat_ind) = [];
    
elseif isfield(cfg,'rmcfg') && strcmpi(cfg.rmcfg,'yes') % Remove hold-overs from old cfg that are too lare to keep around in newer cfgs
    
    varargout{1} = dat;
    
    for iRM = 1:size(dat_ind,2)
        for iSV = 1:numel(cfg.rm_fld)
            varargout{1}.(varargout{1}.data_name{iRM}).cfg.(cfg.rm_fld{iSV}) = [];
            
            rmv_cfg_ind = find(strcmpi(varargout{1}.(varargout{1}.data_name{iRM}).cfg.cfg_lab,cfg.rm_fld{iSV}));
            varargout{1}.(varargout{1}.data_name{iRM}).cfg.cfg_lab(rmv_cfg_ind) = [];
            
        end
    end
       
elseif isfield(fcfg,'save') && strcmpi(fcfg.save,'yes') % Save data 
    
    sve_str = {}; if isfield(dat,'sve_fld'); sve_str = fieldnames(dat.sve_fld); for iSV = 1:numel(sve_str); eval([sve_str{iSV} ' = dat.sve_fld.(sve_str{iSV});']); end; end
    dat.sve_fld = [];
    
    dat_str = dat.data_name; 
    
    cfg_str = {};
    icfg = 1;
    for iDT = 1:numel(dat_str)
        if isfield(dat.(dat_str{iDT}).cfg,'cfg_lab');
            for iCG = 1:numel(dat.(dat_str{iDT}).cfg.cfg_lab); 
                eval([dat_str{iDT} '_'  dat.(dat_str{iDT}).cfg.cfg_lab{iCG} ' = dat.(dat_str{iDT}).cfg.(dat.(dat_str{iDT}).cfg.cfg_lab{iCG});']);
                cfg_str{icfg} = [dat_str{iDT} '_' dat.(dat_str{iDT}).cfg.cfg_lab{iCG}];
                cfg_lab{iDT}{iCG} = dat.(dat_str{iDT}).cfg.cfg_lab{iCG};
                dat.(dat_str{iDT}).cfg.(dat.(dat_str{iDT}).cfg.cfg_lab{iCG}) = [];
                icfg = icfg + 1;
            end
        end
    end
    
    % Figure out Dependencies
    if isfield(cfg,'sve_app'); load(fcfg.filename,'fle_dep'); end
    
    if ~isfield(cfg,'sve_app'); org_fle_dep.str_nme = fcfg.str_nme; elseif isfield(cfg,'sve_app');  org_fle_dep.str_nme =  fle_dep.str_nme; end
    org_fle_dep.sve_fld = sve_str;
    icfg = 1;
    for iDT = 1:numel(dat_str)
        org_fle_dep.(dat_str{iDT}) = {};
        %         if isfield(dat.(dat_str{iDT}).cfg,'cfg_lab'); for iCG = 1:numel(dat.(dat_str{iDT}).cfg.cfg_lab); cfg_lab{iCG} = cfg_str{icfg}; icfg = icfg + 1; end; end
        if isfield(dat.(dat_str{iDT}).cfg,'cfg_lab');
            org_fle_dep.(dat_str{iDT}) = strcat([dat_str{iDT} '_'],cfg_lab{iDT});
            org_fle_dep.(dat_str{iDT})(2,:) = cfg_lab{iDT};            
        end
    end
    
    for iDS = 1:numel(dat_str); eval([dat_str{iDS} ' = dat.(dat_str{iDS});']); dat.(dat_str{iDS})=[]; end
    
    if ~isfield(cfg,'sve_app')
        fle_dep = org_fle_dep;
        if ~isempty(sve_str) && ~isempty(cfg_str); sve_var = [dat_str(:)' sve_str(:)' cfg_str(:)']; elseif ~isempty(cfg_str); sve_var = [dat_str(:)' cfg_str(:)']; elseif ~isempty(sve_str); sve_var = [dat_str(:)' sve_str(:)']; else sve_var = dat_str(:)'; end
        save(fcfg.filename,sve_var{:},'fle_dep','-v7.3')
    elseif isfield(cfg,'sve_app')
        new_sve = setdiff(sve_str,fle_dep.sve_fld);
        fle_dep.sve_fld = [fle_dep.sve_fld new_sve];
        cur_nme = fieldnames(fle_dep); cur_nme = cur_nme(3:end)';
        new_nme = setdiff(dat_str,cur_nme);
        for iUP = 1:numel(new_nme); fle_dep = setfield(fle_dep,new_nme{iUP},{}); end
        cfg_cur_nme = fieldnames(fle_dep); cfg_cur_nme = cfg_cur_nme(3:end);
        for iUP = 1:numel(cfg_cur_nme);
            if isfield(org_fle_dep,cfg_cur_nme{iUP})
                cfg_new_nme{iUP} = org_fle_dep.(cfg_cur_nme{iUP});
                if isfield(fle_dep.(cfg_cur_nme{iUP}),cfg_cur_nme{iUP}); cfg_new_nme{iUP} = setdiff(cfg_new_nme{iUP},fle_dep.(cfg_cur_nme{iUP})(1,:)); end
                if ~isempty(cfg_new_nme{iUP});
                    if ~isempty(fle_dep.(cfg_cur_nme{iUP}))
                        fle_dep.(cfg_cur_nme{iUP}) = [fle_dep.(cfg_cur_nme{iUP})(1,:) cfg_new_nme{iUP}; fle_dep.(cfg_cur_nme{iUP})(2,:) cfg_new_nme{iUP}(2,:)];
                    else
                        fle_dep.(cfg_cur_nme{iUP}) = [cfg_new_nme{iUP}(1,:);cfg_new_nme{iUP}(2,:)];
                    end
                end
            end
        end
        cfg_new_nme = [cfg_new_nme{:}];
        
        if strcmpi(cfg.sve_app,'app_lst') && isfield(cfg,'app_lst') % - EJK need to add
            
        elseif strcmpi(cfg.sve_app,'app_new')
            if ~isempty(cfg_new_nme); sve_var = [new_nme(:)' sve_str(:)' [cfg_new_nme(1,:)]]; else sve_var = [new_nme(:)' sve_str(:)' cfg_new_nme]; end
            save(fcfg.filename,sve_var{:},'fle_dep','-append','-v7.3')
        elseif strcmpi(cfg.sve_app,'app_all')
            sve_var = [dat_str(:)' sve_str(:)' cfg_str(:)'];
            save(fcfg.filename,sve_var{:},'fle_dep','-append','-v7.3')
        end
    end
    
elseif isfield(fcfg,'load') && strcmpi(fcfg.load,'yes') % Load data
    
    if ~isfield(fcfg,'sub_fld'); fcfg.sub_fld = {'all'}; end
    
    load(fcfg.file,'fle_dep');
    varargout{1}           = struct;
    varargout{1}.data_name = {};
    
    sve_nme = fle_dep.sve_fld;
    
    if isfield(fcfg,'sub_fld') && isnumeric(fcfg.sub_fld)
        fld_nme = fieldnames(fle_dep); fld_nme(1:2) = [];
        fcfg.sub_fld = fld_nme(fcfg.sub_fld)';
    end
    
    sub_sve = sum(cell2mat(cellfun(@(x) strcmpi(x,sve_nme),fcfg.sub_fld,'UniformOutput',0)),2);
    for iSV = 1:numel(sve_nme) % sve_fld
        if logical(sub_sve(iSV)) | strcmpi(fcfg.sub_fld,'all')
            sve_tmp                             = load(fcfg.file,sve_nme{iSV});
            varargout{1}.sve_fld.(sve_nme{iSV}) = sve_tmp.(sve_nme{iSV});
        end
    end
    
    % data_fld
    dat_nme = fieldnames(fle_dep);
    dat_sve = sum(cell2mat(cellfun(@(x) strcmpi(x,dat_nme),fcfg.sub_fld,'UniformOutput',0)),2);
    for iDT = 3:numel(dat_nme)
        if logical(dat_sve(iDT)) | strcmpi(fcfg.sub_fld,'all')
            dat_tmp                       = load(fcfg.file,dat_nme{iDT});
            varargout{1}.(dat_nme{iDT})   = dat_tmp.(dat_nme{iDT});
            varargout{1}.data_name{end+1} = dat_nme{iDT};
        end
        
        if ~isempty(fle_dep.(dat_nme{iDT}))
            cfg_nme = fieldnames(fle_dep);
            cfg_sve = sum(cell2mat(cellfun(@(x) strcmpi(x,cfg_nme),fcfg.sub_fld,'UniformOutput',0)),2);
            for iCG = 1:numel(cfg_nme) % cfg
                if  logical(cfg_sve(iCG)) | strcmpi(fcfg.sub_fld,'all')
                    for iCG = 1:numel(fle_dep.(dat_nme{iDT})(1,:))
                        cfg_tmp                                        = load(fcfg.file,fle_dep.(cfg_nme{iDT}){1,iCG});
                        varargout{1}.(dat_nme{iDT}).cfg.(fle_dep.(cfg_nme{iDT}){2,iCG}) = cfg_tmp.(fle_dep.(cfg_nme{iDT}){1,iCG});
                    end
                end
            end
        end
    end
    
else % Normal function
    
    for iRN = 1:size(dat_ind,2) % Loop over fields
        cfg.iRN = iRN;
        
        if isfield(cfg,'pow_phs') || isfield(cfg,'avg_frq'); % if dealing with frequencies & complex numbers
            org_cpl = dat.(dat.data_name{dat_ind(1,cfg.iRN)}).fourierspctrm; rtn_cpl = 1;
            
            if isfield(cfg,'pow_phs')
                if strcmpi(cfg.pow_phs,'pow'); dat.(dat.data_name{dat_ind(1,cfg.iRN)}).fourierspctrm = abs(dat.(dat.data_name{dat_ind(1,cfg.iRN)}).fourierspctrm); elseif strcmpi(cfg.pow_phs,'phs'); dat.(dat.data_name{dat_ind(1,cfg.iRN)}).fourierspctrm = angle(dat.(dat.data_name{dat_ind(1,cfg.iRN)}).fourierspctrm); end
                rtn_for = 0;
            end
            
            if isfield(cfg,'avg_frq')
                rtn_for = 1;
                rtn_trl = 1;
                [tmp,cfg.avg_frq(1)] = min(abs(dat.(dat.data_name{dat_ind(1,cfg.iRN)}).freq-cfg.avg_frq(1)));
                [tmp,cfg.avg_frq(2)] = min(abs(dat.(dat.data_name{dat_ind(1,cfg.iRN)}).freq-cfg.avg_frq(2)));
                dat.(dat.data_name{dat_ind(1,cfg.iRN)}).trial = cellfun(@squeeze,num2cell(squeeze(mean(dat.(dat.data_name{dat_ind(1,cfg.iRN)}).fourierspctrm(:,:,cfg.avg_frq(1):cfg.avg_frq(2),:),3)),[2 3]),'UniformOutput',0);
                dat.(dat.data_name{dat_ind(1,cfg.iRN)})       = rmfield(dat.(dat.data_name{dat_ind(1,cfg.iRN)}),'fourierspctrm');
                dat.(dat.data_name{dat_ind(1,cfg.iRN)}).time  = repmat({dat.(dat.data_name{dat_ind(1,cfg.iRN)}).time},1,numel(dat.(dat.data_name{dat_ind(1,cfg.iRN)}).trialinfo));
            else rtn_trl = 0;
            end
            
        else rtn_cpl = 0; rtn_for = 0; rtn_trl = 0;
        end
        
        if isfield(cfg,'specific'); % To allow fields to be passed to correct data structures
            for iMU = 1:numel(cfg.specific(2,:)); if isnumeric(cfg.specific{2,iMU}); cfg.specific{2,iMU} = num2cell(cfg.specific{2,iMU}); end; end
            for iMU = 1:numel(cfg.specific(1,:)); cfg.(cfg.specific{1,iMU}) = fcfg.(cfg.specific{1,iMU}){cfg.specific{2,iMU}{iRN}}; end;
        end
        
        if isfield(cfg,'multi'); % Set up running multiple values on the same data
            for iMUS = 1:numel(cfg.multi(2,:)); if isnumeric(cfg.multi{2,iMUS}); cfg.multi{2,iMUS} = num2cell(cfg.multi{2,iMUS}); end; end
%             for iMUS = 1:numel(cfg.multi(1,:)); cfg.(cfg.multi{1,iMUS}) = fcfg.(cfg.multi{1,iMUS}){cfg.multi{2,iMUS}{iRN}}; end;
        end
        
        if ~isfield(cfg,'multi'); multi = 1; else multi = numel(cfg.multi{2,1}); end % This is to set up the ability to pass different values to the same function (like band-stop filtering)      
       
        if isfield(fcfg,'dat_sup'); dat_sup = dat.(dat.data_name{dat_ind(1,cfg.iRN)}).cfg.(fcfg.dat_sup); end % Access relevant fields of the cfg to supplemental data
        if isfield(fcfg,'dat_rep');
            rtn_dta = 1;
            if ~isfield(dat.(dat.data_name{dat_ind(1,cfg.iRN)}),'fourierspctrm')
                org_dta = dat.(dat.data_name{dat_ind(1,cfg.iRN)}).trial;
                dat.(dat.data_name{dat_ind(1,cfg.iRN)}).trial = getfield(dat.(dat.data_name{dat_ind(1,cfg.iRN)}).cfg,fcfg.dat_rep{:});
            else
                org_dta = dat.(dat.data_name{dat_ind(1,cfg.iRN)}).fourierspctrm;
                dat.(dat.data_name{dat_ind(1,cfg.iRN)}).fourierspctrm = getfield(dat.(dat.data_name{dat_ind(1,cfg.iRN)}).cfg,fcfg.dat_rep{:});
            end
        else rtn_dta = 0;
        end
            
        if isvar('dat'); dat_tmp = dat; end % In case data needs to be loaded or trials defined
        
        for iMU = 1:multi
            
            if isfield(cfg,'alt_eve'); alt_ret = 1; tmp_org_eve = dat.(dat.data_name{dat_ind(1,cfg.iRN)}).trialinfo; dat.(dat.data_name{dat_ind(1,cfg.iRN)}).trialinfo = dat.(dat.data_name{dat_ind(1,cfg.iRN)}).cfg.alt_eve.(cfg.alt_eve); else alt_ret = 0; end
            if isfield(cfg,'alt_dat'); alt_rmv = {}; for iAL = 1:numel(cfg.alt_dat); alt_rmv = [alt_rmv cfg.alt_dat{iAL}]; dat.(dat.data_name{dat_ind(1,cfg.iRN)}).cfg.alt_dat.(cfg.alt_dat{iAL}) = dat.(dat.data_name{dat_ind(1,cfg.iRN)}).cfg.(cfg.alt_dat{iAL}); end; else alt_rmv = []; end
            
            if isfield(cfg,'multi'); for iMUS = 1:numel(cfg.multi(1,:)); cfg.(cfg.multi{1,iMUS}) = fcfg.(cfg.multi{1,iMUS}){cfg.multi{2,iMUS}{iMU}}; end; end % Set up running multiple values
            
            if nargin==2 % Trl functions
                dat_tmp = feval(func,cfg);
                            
            elseif size(dat_ind,1) == 1 % Typical input (function, cfg, and dataset) for fieldtrip
                
                if strcmpi(char(func),'ft_load_nsx') || (strcmpi(char(func),'ft_preprocessing') && (isfield(cfg,'continuous') || isfield(cfg,'dataset'))) || strcmpi(char(func),'ft_definetrial') % For loading new datafiles
                    dat_tmp = ft_func_v2(func,cfg);
                elseif isfield(cfg,'empty') && strcmpi(cfg.empty,'yes')
                    ft_func_v2(func,cfg,dat.(dat.data_name{dat_ind(1,cfg.iRN)}));
                elseif isfield(fcfg,'dat_sup') 
                    dat_tmp = ft_func_v2(func,cfg,dat_sup,dat.(dat.data_name{dat_ind(1,cfg.iRN)})); % If need to pass additional information from dat cfg
                elseif strcmpi(cfg.data_new,'cfg')
                    dat_tmp = feval(func,cfg,dat.(dat.data_name{dat_ind(1,cfg.iRN)}));
                else
                    dat_tmp = ft_func_v2(func,cfg,dat.(dat.data_name{dat_ind(1,cfg.iRN)})); % For normal processing of single data structures
                end
                
            elseif size(dat_ind,1) == 2 % Added if needed to concatnate structures (must concatenate only two at a time) or other 2 data structures (example: ft_er_pac)
                dat_tmp = ft_func_v2(func,cfg,dat.(dat.data_name{dat_ind(1,cfg.iRN)}),dat.(dat.data_name{dat_ind(2,cfg.iRN)}));
                
            else
                error('Only up to 2 datasets are supported at this time')
            end
            
            if alt_ret; dat_tmp.trialinfo = tmp_org_eve; end
            if rtn_cpl; dat_tmp.fourierspctrm = org_cpl; end
            if rtn_for; dat_tmp = rmfield(dat_tmp,'trial'); end
            if ~isempty(alt_rmv); dat_tmp.cfg = rmfield(dat_tmp.cfg,alt_rmv); end
            if rtn_dta;
                if ~isfield(dat_tmp,'fourierspctrm')
                    dat_tmp.trial = org_dta;
                else
                    dat_tmp.fourierspctrm = org_dta;
                end
            end
            
            % Add temporary data to real data
            if (isempty(cfg.data_new) || strcmpi(cfg.data_new,'no')) && ~isfield(cfg,'empty') % If not adding an additional dataset to structure
                if isfield(dat_tmp,'trial')
                    new_name = findname(func,cfg,dat_ind,dat,dat_tmp);
                    
                    if ~strcmpi(new_name,dat.data_name{dat_ind(1,cfg.iRN)}) % Change name of datafield
                        dat.(new_name) = dat_tmp;
                        dat = rmfield(dat,dat.data_name{dat_ind(1,cfg.iRN)}); % If new_name is the same as the old_name, don't change anything
                    else
                        dat.(new_name) = dat_tmp;
                    end
                    
                    dat.data_name{dat_ind(1,cfg.iRN)} = new_name;
                else
                    new_name = findname(func,cfg,dat_ind,dat,dat_tmp);
                    if ~any(strcmpi(new_name,dat.data_name)); if isfield(dat_tmp,'data_name'); dat.data_name = [dat.data_name new_name]; else dat.data_name = new_name; end; end
                    dat.(new_name) = dat_tmp;
                end
            elseif isfield(fcfg,'data_new') && strcmpi(fcfg.data_new,'cfg') && isfield (fcfg,'cfg_lab')
                new_name = findname(func,cfg,dat_ind,dat,dat_tmp); 
                dat.(new_name).cfg.(fcfg.cfg_lab) = dat_tmp;
                if ~isfield(dat.(new_name).cfg,'cfg_lab'); dat.(new_name).cfg.cfg_lab = {fcfg.cfg_lab}; else dat.(new_name).cfg.cfg_lab = [dat.(new_name).cfg.cfg_lab fcfg.cfg_lab]; end
            elseif strcmpi(cfg.data_new,'yes') && ~isfield(cfg,'empty'); % If adding an additional dataset to structure, figure out proper name
                new_name = findname(func,cfg,dat_ind,dat,dat_tmp); dat.(new_name) = dat_tmp;
                if ~any(strcmpi(dat.data_name,new_name)); dat.data_name = [dat.data_name new_name]; end % Make sure name does not already exist (usually during load)
            end
            
        end
        
        if rtn_trl && (~isfield(cfg,'empty') || ~strcmpi(cfg.empty,'yes')); 
            dat.(new_name).trialinfo = dat.(new_name).cfg.orig_trl(:,4); 
        end
        
    end
    
    varargout{1} = dat;
    
end

end
%% subfunctions
% Very important name control of data subfields
function new_name = findname(nfunc,ncfg,ndat_ind,ndat,ndat_tmp) %% - EJK NEED TO ADD ABILITY TO NAME ON NEW STRUCTURES

org_name = ndat.data_name{ndat_ind(1,ncfg.iRN)};

if isfield(ncfg,'new_suffix') && iscell(ncfg.new_suffix)
    new_suf = ncfg.new_suffix{ncfg.iRN};
elseif isfield(ncfg,'new_suffix') && ~iscell(ncfg.new_suffix)
    new_suf = ncfg.new_suffix;
elseif isfield(ncfg,'new_name') && iscell(ncfg.new_name)
    new_nme = ncfg.new_name{ncfg.iRN};
elseif isfield(ncfg,'new_name') && ~iscell(ncfg.new_name)
    new_nme = ncfg.new_name;
end

if strcmpi(char(nfunc),'ft_preprocessing') && isfield(ncfg,'continuous') % This names the new datafields when loading new data
    new_name = mmil_spec_char(ncfg.data_name{ndat_ind(1,ncfg.iRN)},{'-','.'});
elseif isfield(ncfg,'new_suffix')
    new_name = [org_name '_' new_suf];
elseif isfield(ncfg,'new_name')
    new_name = new_nme;
else
    new_name = org_name;
end

end