%% FT_REDEFINE_EVENTS
% This function will change the events into new events (combining events
% into larger events). It saves the original list of events in the data
% structure, and will return that list if the option is chosen.
% 
% cfg.return_events = 0 if combining events, 1 if returning original event list
% cfg.old_events = cell array of vectors of events to combine
% cfg.new_events = cell array of numbers matching the vectors above
% cfg.crt_alt_eve = create field of alternate events in cfg
% 
% TO BE ADDED: How to move the same event into 2 new events, struct of
% event masks
% 
% Created by Erik Kaestner (3-13-14) ekaestne@ucsd.edu

function data = ejk_redefine_events(cfg,data)

if (~isfield(cfg,'return_events')) 
    error('Need to specify whether to return the original data or create new data')
elseif cfg.return_events == 1 && (~isfield(data.cfg,'orig_trialinfo')) 
    error('Cannot return original events, ft_redefine_events has not been used on this data structure')
end

if cfg.return_events == 0;
    
    if isfield(cfg,'new_events')
        
        if isfield(cfg,'crt_alt_eve')
            
            if isfield(cfg,'use_alt_eve')
                red_ind = cell(1,length(cfg.old_events));
                for ired = 1:length(cfg.old_events)
                    red_ind{ired} = ismember(data.cfg.alt_eve.(cfg.use_alt_eve),cfg.old_events{ired});
                end
                
                data.cfg.alt_eve.orig_trialinfo.(cfg.use_alt_eve) = data.trialinfo;
                data.cfg.alt_eve.(cfg.crt_alt_eve) = zeros(numel(data.cfg.alt_eve.(cfg.use_alt_eve)),1);
                for ired = 1:length(cfg.new_events)
                    data.cfg.alt_eve.(cfg.crt_alt_eve)(red_ind{ired}) = cfg.new_events{ired};
                end
            else
                red_ind = cell(1,length(cfg.old_events));
                for ired = 1:length(cfg.old_events)
                    red_ind{ired} = ismember(data.trialinfo,cfg.old_events{ired});
                end
                
                data.cfg.orig_trialinfo = data.trialinfo;
                data.cfg.alt_eve.(cfg.crt_alt_eve) = zeros(numel(data.trialinfo),1);
                for ired = 1:length(cfg.new_events)
                    data.cfg.alt_eve.(cfg.crt_alt_eve)(red_ind{ired}) = cfg.new_events{ired}';
                end
            end
        elseif isfield(cfg,'use_alt_eve')
            red_ind = cell(1,length(cfg.old_events));
            for ired = 1:length(cfg.old_events)
                red_ind{ired} = ismember(data.cfg.alt_eve.(cfg.use_alt_eve),cfg.old_events{ired});
            end
            
            data.cfg.alt_eve.orig_trialinfo.(cfg.use_alt_eve) = data.trialinfo;
            data.cfg.alt_eve.(cfg.use_alt_eve) = zeros(numel(data.cfg.alt_eve.(cfg.use_alt_eve)),1);
            for ired = 1:length(cfg.new_events)
                data.cfg.alt_eve.(cfg.use_alt_eve)(red_ind{ired}) = cfg.new_events{ired};
            end
        else
            if isfield(data.cfg,'orig_trialinfo') && ~isempty(data.cfg.orig_trialinfo) && ~isfield(cfg,'crt_alt_eve') && ~isfield(cfg,'use_alt_eve') && ~isfield(cfg,'plt_hck')
                error('Events already renamed. Either return old events before proceeding or use crt_alt_eve or use_alt_eve');
                %             red_ind = cell(1,length(cfg.old_events));
                %             for ired = 1:length(cfg.old_events)
                %                 red_ind{ired} = ismember(data.cfg.orig_trialinfo,cfg.old_events{ired});
                %             end
                %
                %             for ired = 1:length(cfg.new_events)
                %                 data.trialinfo(red_ind{ired}) = cfg.new_events{ired};
                %             end
                
            else
                
                red_ind = cell(1,length(cfg.old_events));
                for ired = 1:length(cfg.old_events)
                    red_ind{ired} = ismember(data.trialinfo,cfg.old_events{ired});
                end
                
                data.cfg.orig_trialinfo = data.trialinfo;
                
                data.trialinfo = zeros(numel(data.trialinfo),1);
                for ired = 1:length(cfg.new_events)
                    data.trialinfo(red_ind{ired}) = cfg.new_events{ired};
                end
                
            end
        end
        
    elseif isfield(cfg,'eve_rul')
        
        if isfield(cfg,'fle_nme') && ~isempty(string_find({cfg.fle_nme{1,1}},'.csv'))
            
        elseif isfield(cfg,'fle_nme') && ~isempty(string_find({cfg.fle_nme{1,1}},'.mat'))
            ttt = load(cfg.fle_nme{1,1}); ttt_fld = fieldnames(ttt); ttt = ttt.(ttt_fld{1});
            if size(cfg.fle_nme,1)==4; ttt = cat(1,ttt{cfg.fle_nme{4,1}}); end
            ttt = ttt(cfg.fle_nme{2,1}:end,cfg.fle_nme{3,1}:end);
        elseif isfield(cfg,'fld_nme')
            ttt = data.cfg.alt_eve.(cfg.fld_nme);
        end
        
        
        if strcmpi(cfg.eve_rul,'med_spl')
            
        elseif strcmpi(cfg.eve_rul,'eve_sql')
            eve_hld_one = ttt(data.cfg.alt_eve.(cfg.use_alt_eve)==cfg.old_events{1}(1));
            eve_hld_two = ttt(data.cfg.alt_eve.(cfg.use_alt_eve)==cfg.old_events{1}(2));
            
            eve_num_one = find(data.cfg.alt_eve.(cfg.use_alt_eve)==cfg.old_events{1}(1));
            eve_num_two = find(data.cfg.alt_eve.(cfg.use_alt_eve)==cfg.old_events{1}(2));
            
            rmv_ind_one = [];
            rmv_ind_two = [];
            
            [~,p_hld] = ttest2(eve_hld_one,eve_hld_two);
            while p_hld < 0.1
                
                if mean(eve_hld_one) > mean(eve_hld_two)
                    [~,max_ind] = max(eve_hld_one); [~,min_ind] = min(eve_hld_two);
                    
                    rmv_ind_one = [rmv_ind_one eve_num_one(max_ind)]; rmv_ind_two = [rmv_ind_two eve_num_two(min_ind)];
                    
                    eve_hld_one(max_ind) = []; eve_num_one(max_ind) = [];
                    eve_hld_two(min_ind) = []; eve_num_two(min_ind) = [];
                    
                elseif mean(eve_hld_one) < mean(eve_hld_two)
                    [~,max_ind] = max(eve_hld_two); [~,min_ind] = min(eve_hld_one);
                    
                    rmv_ind_one = [rmv_ind_one eve_num_one(min_ind)]; rmv_ind_two = [rmv_ind_two eve_num_two(max_ind)];
                    
                    eve_hld_one(min_ind) = []; eve_num_one(min_ind) = [];
                    eve_hld_two(max_ind) = []; eve_num_two(max_ind) = [];
                    
                end
                
                [~,p_hld] = ttest2(eve_hld_one,eve_hld_two);
            end
            
            data.cfg.alt_eve.(cfg.use_alt_eve)(rmv_ind_one) = 0;
            data.cfg.alt_eve.(cfg.use_alt_eve)(rmv_ind_two) = 0;
            
        else
            
        end
        
    end
        
else
    
    if isfield(data.cfg,'orig_trialinfo') && ~isempty(data.cfg.orig_trialinfo)
        data.trialinfo = data.cfg.orig_trialinfo;
    else
        error('No changed events to move back to the original')
    end
    
end

end