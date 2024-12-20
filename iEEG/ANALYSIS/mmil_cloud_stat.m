function rtn_dat = mmil_cloud_stat(fcfg,dat)

if ~isfield(fcfg,'loc'); fcfg.loc = 'cluster'; end
if ~isfield(fcfg,'chn'); fcfg.chn = 1:numel(dat.label); end

rtn_dat = dat;

if strcmpi(fcfg.loc,'cluster') || strcmpi(fcfg.loc,'local')
    
    bsc_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/tmp_stat';
    
    for iTT = 1:numel(fcfg.stt_fnc)
        
        %% Find old original stats
        if isfield(rtn_dat.cfg,'alt_stt')
            old_stt = fieldnames(rtn_dat.cfg.alt_stt);
        end
        
        if strcmpi(fcfg.loc,'cluster')
            
            %% create holding folder
            fld = 0;
            nme_cnt = 1;
            while ~fld
                if ~exist([bsc_fld '/' fcfg.fld_nme '_' fcfg.stt_fnc{iTT} '_' num2str(nme_cnt)],'dir')
                    fld = 1;
                else
                    nme_cnt = nme_cnt+1;
                end
            end
            mkdir([bsc_fld '/' fcfg.fld_nme '_' fcfg.stt_fnc{iTT} '_' num2str(nme_cnt) '/' 'hold']);
            mkdir([bsc_fld '/' fcfg.fld_nme '_' fcfg.stt_fnc{iTT} '_' num2str(nme_cnt) '/' 'touch']);
            mkdir([bsc_fld '/' fcfg.fld_nme '_' fcfg.stt_fnc{iTT} '_' num2str(nme_cnt) '/' 'product']);
            
            %% save data to holding folder
            parfor iCH = fcfg.chn
                
                stt_dat = [];
                
                tmp_hld = [fcfg.fld_nme '_chn_' mmil_spec_char(dat.label{iCH},{'-' ' ' '>' '(' ')'})];
                if numel(tmp_hld) > 54;
                    tmp_hld = [tmp_hld(1:49) '_' num2str(iCH)];
                end
                stt_dat.data_name{1} = tmp_hld;
                stt_dat.(stt_dat.data_name{1}) = dat;
                
                cfg = [];
                cfg.channel = dat.label{iCH};
                stt_dat = ft_func(@ft_preprocessing,cfg,stt_dat);
                
                parsave([bsc_fld '/' fcfg.fld_nme '_' fcfg.stt_fnc{iTT} '_' num2str(nme_cnt) '/' 'hold' '/' stt_dat.data_name{1} '.mat'],stt_dat)
                
            end          
            
        end
        
        %% make scripts
        if strcmpi(fcfg.loc,'cluster')
            
            scp_fld  = '/home/ekaestne/batchdirs';
            batchdir = [scp_fld '/'  fcfg.fld_nme '_' fcfg.stt_fnc{iTT} '_' num2str(nme_cnt)];
            mkdir(batchdir);
            
            template_file = which(fcfg.stt_fnc{iTT});
            
            for iST = fcfg.chn
                % create scriptlist
                tmp_hld = [fcfg.fld_nme '_chn_' mmil_spec_char(dat.label{iST},{'-' ' ' '>' '(' ')'})];
                if numel(tmp_hld) > 54;
                    tmp_hld = [tmp_hld(1:49) '_' num2str(iST)];
                end
                
                if ~exist(sprintf('%s/scriptlist.txt',batchdir),'file')
                    fid =  fopen(sprintf('%s/scriptlist.txt',batchdir),'w');
                    fprintf(fid,'job_%i_%s',iST,tmp_hld);
                else
                    fid =  fopen(sprintf('%s/scriptlist.txt',batchdir),'a');
                    fprintf(fid,'\njob_%i_%s',iST,tmp_hld);
                end
                fclose(fid);
                
                % creates script
                fid     = fopen(template_file);
                temp_m  = fscanf(fid,'%c',Inf);
                temp_m   = strrep(temp_m, '%', '%%');
                temp_m   = strrep(temp_m, '\', '\\');
                
                str_m1 = ['load(''' bsc_fld '/' fcfg.fld_nme '_' fcfg.stt_fnc{iTT} '_' num2str(nme_cnt) '/' 'hold' '/' tmp_hld '.mat'');'];
                str_m2 = ['stt_dat = stt_dat.(stt_dat.data_name{1});'];
                str_m3 = ['addpath(genpath(''/space/mdeh4/1/halgdev/projects/mmilanguage/matlab''));'];
                end_m1 = ['stt_dat = stt_dat.cfg.alt_stt;'];
                end_m2 = ['save(''' bsc_fld '/' fcfg.fld_nme '_' fcfg.stt_fnc{iTT} '_' num2str(nme_cnt) '/' 'product' '/' tmp_hld '.mat'',''stt_dat'');'];
                add_m  = ['system(''touch ' bsc_fld '/' fcfg.fld_nme '_' fcfg.stt_fnc{iTT} '_' num2str(nme_cnt) '/' 'touch' '/' tmp_hld ''');'];
                
                new_m   = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s',str_m1,str_m2,str_m3,temp_m,end_m1,end_m2,add_m);
                
                fclose(fid);
                if template_file(regexp(template_file,'\.')+1:end) == 'm';
                    fid =fopen(sprintf('%s/job_%i_%s.m',batchdir,iST,tmp_hld),'w');
                else
                    fid =fopen(sprintf('%s/job_%i_%s.csh',batchdir,iST,tmp_hld),'w');
                end
                fprintf(fid,char(new_m(:)));
                fclose(fid);
            end
            
        elseif strcmpi(fcfg.loc,'local')
            
            % Get script
            template_file = which(fcfg.stt_fnc{iTT});
            
            fid     = fopen(template_file);
            temp_m  = fscanf(fid,'%c',Inf);
            temp_m   = strrep(temp_m, '%', '%%');
            temp_m   = strrep(temp_m, '\', '\\');
            
            add_stt_hld_one = strfind(temp_m,'cfg.add_stt');
            add_stt_hld_two = strfind(temp_m(add_stt_hld_one(1):end),'''');
            add_stt_hld = temp_m(add_stt_hld_two(1)+add_stt_hld_one:add_stt_hld_two(2)+add_stt_hld_one-2);
            
            % Add in new script
            str_a1 = ['function dat = tmp_' fcfg.stt_fnc{iTT} '(dat)'];
            
            str_a2   = ['dat.cfg.alt_stt.(''' add_stt_hld ''').prob  = nan(numel(dat.label),numel(dat.time{1}));'];
            str_a2_5 = ['prb_stt_hld = nan(' num2str(numel(fcfg.chn)) ',numel(dat.time{1}));'];
            
            str_fnd_beg = strfind(temp_m,'cfg.alt_eve ');
            str_fnd_end = strfind(temp_m,';');
            str_fnd_end = str_fnd_end(str_fnd_end>str_fnd_beg);
            str_fnd_end = str_fnd_end(dsearchn(str_fnd_end',str_fnd_beg));
            
            str_a3 = 'pcfg = [];';
            str_a4 = ['p' temp_m(str_fnd_beg:str_fnd_end)];
            
            str_a5 = '[typ,dta,lbl,eve_hld,tme] = mmil_prep_stat(pcfg,dat);';
            str_a5_5 = [ 'chn_use = ' '[' num2str(fcfg.chn) '];' ];  
            
            str_a6 = 'parfor iC = 1:numel(chn_use)'; %['parfor iC = [' num2str(fcfg.chn) ']']; % fcfg.chn ['parfor iC = 1:numel(dat.label)'];

            str_a7 = 'if ~(sum(isnan(dta(chn_use(iC),1,:)))==numel(squeeze(dta(chn_use(iC),1,:))))';
            
            str_a8  = 'cfg = [];';
            str_a9  = 'cfg.typ = typ;';
            str_a10  = 'cfg.dta = dta(chn_use(iC),:,:);';
            str_a11 = 'cfg.lbl = lbl(chn_use(iC));';
            str_a12 = 'cfg.eve_hld = eve_hld;';
            str_a13 = 'cfg.tme = tme;';
                        
            temp_m  = strrep(temp_m,'stt_dat','prb_stt_hld(iC,:)');
            
%             str_e1 = ['prb_stt_hld(iC,:) = chn_dat{iC}.cfg.alt_stt.' add_stt_hld '.prob;'];
            str_e1 = 'end';
            str_e2 = 'end';
            str_e2_5 = 'for iC = 1:numel(chn_use)';
            str_e3 = ['dat.cfg.alt_stt.(''' add_stt_hld ''').prob(chn_use(iC),:)  = prb_stt_hld(iC,:);'];
            str_e3_5 = 'end';
            str_e4 = ['dat.cfg.alt_stt.(''' add_stt_hld ''').label = dat.label;'];
            str_e5 = ['dat.cfg.alt_stt.(''' add_stt_hld ''').time  = dat.time;'];
            
            
            new_m   = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n', ...
                str_a1,str_a2,str_a2_5,str_a3,str_a4,str_a5,str_a5_5,str_a6,str_a7,str_a8,str_a9,str_a10,str_a11,str_a12,str_a13,temp_m,str_e1,str_e2,str_e2_5,str_e3,str_e3_5,str_e4,str_e5);
           
            fclose(fid);
            if numel([add_stt_hld '_' fcfg.fld_nme '_func.m'])-2>63; sub = numel([add_stt_hld '_' fcfg.fld_nme '_func.m'])-2-62; else sub = 0; end            
            fid =fopen(['/space/mdeh4/1/halgdev/projects/mmilanguage/matlab/scratch/stt_hld/' add_stt_hld '_' fcfg.fld_nme(1:end-sub) '_func.m' ],'w');
            fprintf(fid,char(new_m(:)));
            fclose(fid);
            
        end
        
        %% Submit Jobs & wait
        if strcmpi(fcfg.loc,'local')
            
            addpath('/space/mdeh4/1/halgdev/projects/mmilanguage/matlab/scratch/stt_hld')
            eval(['rtn_dat = ' add_stt_hld '_'  fcfg.fld_nme(1:end-sub) '_func(rtn_dat);' ])
            
        elseif strcmpi(fcfg.loc,'cluster')
            
            fid =fopen([batchdir '/' 'run_cluster.csh'],'w');
            fprintf(fid,'#!/bin/csh -f\n');
            fprintf(fid,char(['ssh ekaestne@mmilcluster ''qmatjobs ' fcfg.fld_nme '_' fcfg.stt_fnc{iTT} '_' num2str(nme_cnt) '''']));
            fclose(fid);
                  
            system(['sh ' batchdir '/' 'run_cluster.csh'])
            
            % Wait till all stats run
            wait = 1;
            while wait
                num_dat = numel(dir([bsc_fld '/' fcfg.fld_nme '_' fcfg.stt_fnc{iTT} '_' num2str(nme_cnt) '/' 'hold' '/']));
                num_fin = numel(dir([bsc_fld '/' fcfg.fld_nme '_' fcfg.stt_fnc{iTT} '_' num2str(nme_cnt) '/' 'touch' '/']));
                if num_dat == num_fin; wait = 0; end
                pause(10)
            end
            
        end
                
        %% Find new stats and put them together into rtn_dat
        if strcmpi(fcfg.loc,'cluster')
            
            tmp_hld = [fcfg.fld_nme '_chn_' mmil_spec_char(dat.label{1},{'-' ' ' '>' '(' ')'})];
            if numel(tmp_hld) > 54;
                tmp_hld = [tmp_hld(1:49) '_' num2str(1)];
            end
            
            tmp_stt_dat = load([bsc_fld '/' fcfg.fld_nme '_' fcfg.stt_fnc{iTT} '_' num2str(nme_cnt) '/' 'product' '/' tmp_hld '.mat']);
            new_stt = fieldnames(tmp_stt_dat.stt_dat);
            if isfield(rtn_dat.cfg,'alt_stt')
                rmv_ind = [];
                old_stt = fieldnames(rtn_dat.cfg.alt_stt);
                for iNS = 1:numel(new_stt)
                    stt_ext_ind = find(strcmpi(new_stt{iNS},old_stt));
                    if ~isempty(stt_ext_ind)
                        if isfield(rtn_dat.cfg.alt_stt.(new_stt{stt_ext_ind}),'cfg'); old_num = numel(fieldnames(rtn_dat.cfg.alt_stt.(new_stt{stt_ext_ind}).cfg.script)); else old_num = 0; end
                        if isfield(tmp_stt_dat.stt_dat.(new_stt{iNS}),'cfg') && ~(numel(fieldnames(tmp_stt_dat.stt_dat.(new_stt{iNS}).cfg.script)) > old_num)
                            rmv_ind = [rmv_ind iNS];
                        elseif ~isfield(tmp_stt_dat.stt_dat.(new_stt{iNS}),'cfg') && ~isfield(rtn_dat.cfg.alt_stt.(new_stt{iNS}),'cfg')
                            rmv_ind = [rmv_ind iNS];
                        end
                    end
                end
            end
            if isvar('rmv_ind') && ~isempty(rmv_ind); new_stt(rmv_ind) = []; end
            
            stt_fle = dir([bsc_fld '/' fcfg.fld_nme '_' fcfg.stt_fnc{iTT} '_' num2str(nme_cnt) '/' 'product' '/' '*.mat' ]); stt_fle = {stt_fle(:).name};
            prb_hld = zeros(numel(new_stt),numel(tmp_stt_dat.stt_dat.(new_stt{1}).prob),numel(stt_fle));
            prb_sze = 1:size(prb_hld,1);
            parfor iL = 1:numel(stt_fle);
                
                stt_dat = load([bsc_fld '/' fcfg.fld_nme '_' fcfg.stt_fnc{iTT} '_' num2str(nme_cnt) '/' 'product' '/' stt_fle{iL}]);
                
                label{iL} = stt_dat.stt_dat.(new_stt{1}).label{1};
                
                spc_prb_hld = zeros(numel(new_stt),numel(tmp_stt_dat.stt_dat.(new_stt{1}).prob));
                
                for iST = prb_sze
                    spc_prb_hld(prb_sze(iST),:) = stt_dat.stt_dat.(new_stt{iST}).prob(1,:);
                end
                
                prb_hld(:,:,iL) = spc_prb_hld;
                
            end
            
            for iST = prb_sze
                rtn_dat.cfg.alt_stt.(new_stt{iST}).label = label;
                rtn_dat.cfg.alt_stt.(new_stt{iST}).prob  = squeeze(prb_hld(iST,:,:))';
                rtn_dat.cfg.alt_stt.(new_stt{iST}).time  = tmp_stt_dat.stt_dat.(new_stt{iST}).time;
                rtn_dat.cfg.alt_stt.(new_stt{iST}).cfg   = tmp_stt_dat.stt_dat.(new_stt{iST}).cfg;
            end
            
            system(['rm -r ' bsc_fld '/' fcfg.fld_nme '_' fcfg.stt_fnc{iTT} '_' num2str(nme_cnt)])
            
        end
    end
    
elseif strcmpi(fcfg.loc,'not')
    
    for iTT = 1:numel(fcfg.stt_fnc)
        
        template_file = which(fcfg.stt_fnc{iTT});
        
        spl_pnt = strfind(template_file,'/');
        tmp_fle = template_file(spl_pnt(end)+1:end);
        
        % creates script
        fid     = fopen(template_file);
        temp_m  = fscanf(fid,'%c',Inf);
        temp_m   = strrep(temp_m, '%', '%%');
        temp_m   = strrep(temp_m, '\', '\\');
        
        str_m1 = ['function stt_dat = ' tmp_fle(1:end-2) '_' fcfg.fld_nme '(stt_dat)'];
        add_m  = ['end'];
            
        new_m   = sprintf('%s\n%s\n%s\n%s\n%s',str_m1,temp_m,add_m);
        
        fclose(fid);
        fid =fopen(sprintf('/space/mdeh4/1/halgdev/projects/mmilanguage/tmp_stat/loc_tmp_stt/%s_%s.m',tmp_fle(1:end-2),fcfg.fld_nme),'w');
        fprintf(fid,char(new_m(:)));
        fclose(fid);
        addpath('/space/mdeh4/1/halgdev/projects/mmilanguage/tmp_stat/loc_tmp_stt/')
        
        eval(['rtn_dat = ' tmp_fle(1:end-2) '_' fcfg.fld_nme '(rtn_dat)']);
        
    end
    
end

end

%% OLD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             fle = mmil_readtext([batchdir '/' 'scriptlist.txt']);
%             
%             parfor iRN = 1:numel(fle)
%                 parrun([batchdir '/' fle{iRN} '.m'])
%             end