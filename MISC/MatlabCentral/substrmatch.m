function idx = substrmatch(word,cellarray)
    idx = ~cellfun(@isempty,strfind(word,cellarray));