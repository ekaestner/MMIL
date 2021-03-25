function string = mmil_spec_char(string,character,varargin)

for iCH = 1:numel(character)
    
    if isempty(varargin)
        if ~strcmpi(character{iCH},'_'); rep_chr = {'_'}; else rep_chr = {''}; end
    else
        rep_chr = varargin{:};
    end
    
    spc_ind = strfind(string,character{iCH});
    for iRP = 1:length(spc_ind)
        if ~isempty(spc_ind(iRP))
        string  = [string(1:spc_ind(iRP)-1) rep_chr{1} string(spc_ind(iRP)+1:end)];
        end
    end
end

if ~isempty(string) && strcmpi(string(1),'_'); string(1) = 'x'; end

end