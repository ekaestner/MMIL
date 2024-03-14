function string = mmil_spec_char(string,character,varargin)

if ~iscell(string); string = {string}; sng_str = 1; else sng_str = 0; end

for iST = 1:numel(string)
    for iCH = 1:numel(character)
        
        if isempty(varargin)
            if ~strcmpi(character{iCH},'_'); rep_chr = {'_'}; else rep_chr = {''}; end
        else
            rep_chr = varargin{:};
        end
        
        spc_ind = strfind(string{iST},character{iCH});
        for iRP = 1:length(spc_ind)
            if ~isempty(spc_ind(iRP))
                string{iST}  = [string{iST}(1:spc_ind(iRP)-1) rep_chr{1} string{iST}(spc_ind(iRP)+1:end)];
            end
        end
    end
    
    if ~isempty(string{iST}) && strcmpi(string{iST}(1),'_'); string{iST}(1) = 'x'; end
end

if sng_str; string = string{1}; end

end