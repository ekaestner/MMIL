
function ejk_chk_dir(dir_nme)

if ~exist(dir_nme,'dir')
    mkdir(dir_nme)
end

end