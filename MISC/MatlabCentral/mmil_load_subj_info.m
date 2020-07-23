function ele_txt = mmil_load_subj_info(in__fle,clr_fld,varargin)

ele_txt = mmil_readtext(in__fle,'[]');

if ~isempty(varargin) && strcmpi(varargin{1},'update')
    clr_fle = order_sbj_inf(ele_txt);
    cell2csv(in__fle,clr_fle)
else
    fld_nme = cellfun(@(x) strfind(x,'|'),ele_txt,'uni',0);
    fld_nme = cellfun(@(x,y) x(1:y-2),ele_txt,fld_nme,'uni',0);
    fld_nme(cellfun(@isempty,fld_nme)) = {''};
    
    if isempty(varargin)
        
        ele_txt = ele_txt(string_find(fld_nme,clr_fld));
        ele_txt = ele_txt{1}(strfind(ele_txt{1},' | ')+2:end);
        if strcmpi(ele_txt(1),' '); ele_txt = {ele_txt(2:end)}; end
        
        if ~isempty(strfind(ele_txt{1},';'));
            ele_txt = strsplit(ele_txt{1},';');
            for iC = 1:numel(ele_txt)
                ele_txt{iC} = strsplit(ele_txt{iC},',');
            end
            for iC = 1:numel(ele_txt); for iN = 1:numel(ele_txt{iC}); if strcmpi(ele_txt{iC}{iN},'[]') || strcmpi(ele_txt{iC}{iN},'{}'); ele_txt{iC}{iN} = []; end; end; end;
            for iC = 1:numel(ele_txt);
                for iN = 1:numel(ele_txt{iC});
                    if strcmpi(ele_txt{iC}{iN},'[]') || strcmpi(ele_txt{iC}{iN},'{}');
                        ele_txt{iC}{iN} = [];
                    else
                        if ~isempty(strfind(ele_txt{iC}{iN},'}')) & ~isempty(strfind(ele_txt{iC}{iN},'{'))
                            ele_txt{iC}{iN}(strfind(ele_txt{iC}{iN},'}')) = [];
                            ele_txt{iC}{iN}(strfind(ele_txt{iC}{iN},'{')) = [];
                        elseif ~isempty(strfind(ele_txt{iC}{iN},'}'))
                            ele_txt{iC}{iN}(strfind(ele_txt{iC}{iN},'}')) = [];
                        elseif ~isempty(strfind(ele_txt{iC}{iN},'{'))
                            ele_txt{iC}{iN}(strfind(ele_txt{iC}{iN},'{')) = [];
                        end;
                        if ~isempty(ele_txt{iC}{iN}); tmp = str2num(ele_txt{iC}{iN}); if ~isempty(tmp) & ~isnan(tmp); ele_txt{iC}{iN} = str2num(ele_txt{iC}{iN}); end; end
                    end
                end;
            end;
        
        else
            if ~isempty(strfind(ele_txt{1},','));
                ele_txt = strsplit(ele_txt{1},',');
            end
            for iC = 1:numel(ele_txt);
                if strcmpi(ele_txt{iC},'[]') || strcmpi(ele_txt{iC},'{}'); ele_txt{iC} = [];
                else
                    if ~isempty(strfind(ele_txt{iC},'}')) & ~isempty(strfind(ele_txt{iC},'{'))
                        ele_txt{iC}{iN}(strfind(ele_txt{iC},'}')) = [];
                        ele_txt{iC}{iN}(strfind(ele_txt{iC},'{')) = [];
                    elseif ~isempty(strfind(ele_txt{iC},'}'))
                        ele_txt{iC}(strfind(ele_txt{iC},'}')) = [];
                    elseif ~isempty(strfind(ele_txt{iC},'{'))
                        ele_txt{iC}(strfind(ele_txt{iC},'{')) = [];
                    end;
                    if ~exist([ele_txt{iC} '.m'],'file'); tmp = str2num(ele_txt{iC}); if ~isempty(tmp) & ~isnan(tmp); ele_txt{iC} = str2num(ele_txt{iC}); end; end
                end;
            end
        end
    else
        if strcmpi(varargin{1},'asstring')
            ele_txt = ele_txt(string_find(ele_txt,clr_fld));
            ele_txt = ele_txt{1}(strfind(ele_txt{1},' | ')+2:end);
            if strcmpi(ele_txt(1),' '); ele_txt = {ele_txt(2:end)}; end
        else
            if numel(varargin)==1; varargin{2} = 1; end
            ele_rep_loc = string_find(ele_txt,clr_fld);
            ele_txt{ele_rep_loc}(strfind(ele_txt{ele_rep_loc},' | ')+2:end) = [];
            
            if varargin{2} == 1
                ele_txt{ele_rep_loc} = [ele_txt{ele_rep_loc} ' ' varargin{1}];
            elseif varargin{2} == 0
                ele_txt{ele_rep_loc} = [clr_fld ' | ' varargin{1}];
            end
            
            cell2csv(in__fle,ele_txt)
        end
    end
end

end
