function idx = substrmatch(word,cellarray)
idx = ~cellfun(@(x) strcmp(word, x),cellarray);