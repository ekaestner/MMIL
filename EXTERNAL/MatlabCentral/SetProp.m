function obj=SetProp(obj, prop, std_val, varargin)

S   = struct('type','.','subs',regexp(prop,'\.','split'));
obj = subsasgn(obj, S, std_val);

end