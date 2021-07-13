function pth_new = ejk_fix_path( pth_old )

spc_fnd = strfind( pth_old, ' ');

for iSF = 1:numel(spc_fnd)
    pth_old = [ pth_old(1:spc_fnd(iSF)-1) '\' pth_old(spc_fnd(iSF):end) ];
    spc_fnd = strfind( pth_old, ' ');
end

pth_new = pth_old;

end