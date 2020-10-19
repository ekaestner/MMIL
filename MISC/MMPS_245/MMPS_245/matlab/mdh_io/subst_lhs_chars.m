function lhs = subst_lhs_chars(lhs)

lhs = deblank(lhs);
lhs(find(lhs=='.')) = 'd';
lhs(find(lhs=='[')) = 'l';
lhs(find(lhs==']')) = 'r';

