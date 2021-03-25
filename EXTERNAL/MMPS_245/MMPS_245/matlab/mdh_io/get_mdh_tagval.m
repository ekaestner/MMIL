function tagval = get_mdh_field(mdhStruct,tagname)

tagname = subst_lhs_chars(tagname);

if isfield(mdhStruct,tagname)
  tagval = getfield(mdhStruct,tagname);
else
  tagval = 0;
end