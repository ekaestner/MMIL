plt_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/initial_plot';
cpy_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/initial_plot/repetition_effects';

rmdir(cpy_dir,'s')
mkdir(cpy_dir)

% Find all files with 'repetition' in lfp visual folders
tmp_cpy = strsplit(ls([plt_dir '/' '*' '/' 'Visual' '/' '*lfp*' '/' '*_repetition_differences_*']),'.png');
mkdir([cpy_dir '/' 'vis_lfp'])
cellfun(@(x) copyfile([x '.png'],[cpy_dir '/' 'vis_lfp']),tmp_cpy);
   
% Find all files with 'repetition' in hgp visual folders
tmp_cpy = strsplit(ls([plt_dir '/' '*' '/' 'Visual' '/' '*hgp*' '/' '*_repetition_differences_*']),'.png');
mkdir([cpy_dir '/' 'vis_hgp'])
cellfun(@(x) copyfile([x '.png'],[cpy_dir '/' 'vis_hgp']),tmp_cpy);

% Find all files with 'repetition' in lfp auditory folders
tmp_cpy = strsplit(ls([plt_dir '/' '*' '/' 'Auditory' '/' '*lfp*' '/' '*_repetition_differences_*']),'.png');
mkdir([cpy_dir '/' 'aud_lfp'])
cellfun(@(x) copyfile([x '.png'],[cpy_dir '/' 'aud_lfp']),tmp_cpy);

% Find all files with 'repetition' in lfp auditory folders
tmp_cpy = strsplit(ls([plt_dir '/' '*' '/' 'Auditory' '/' '*hgp*' '/' '*_repetition_differences_*']),'.png');
mkdir([cpy_dir '/' 'aud_hgp'])
cellfun(@(x) copyfile([x '.png'],[cpy_dir '/' 'aud_hgp']),tmp_cpy);

% Find all files with 'repetition' in lfp visual/auditory folders
tmp_cpy = strsplit(ls([plt_dir '/' '*' '/' 'Auditory_Visual' '/' '*lfp*' '/' '*_repetition_differences_*']),'.png');
mkdir([cpy_dir '/' 'aud_vis_lfp'])
cellfun(@(x) copyfile([x '.png'],[cpy_dir '/' 'aud_vis_lfp']),tmp_cpy);

% Find all files with 'repetition' in lfp visual/auditory folders
tmp_cpy = strsplit(ls([plt_dir '/' '*' '/' 'Auditory_Visual' '/' '*hgp*' '/' '*_repetition_differences_*']),'.png');
mkdir([cpy_dir '/' 'aud_vis_hgp'])
cellfun(@(x) copyfile([x '.png'],[cpy_dir '/' 'aud_vis_hgp']),tmp_cpy);

% Find all location files with 'repetition' in lfp
tmp_cpy = strsplit(ls([plt_dir '/' '*' '/' 'initial_location' '/' '*lfp*' '/' '*_rep_*']),'.png');
mkdir([cpy_dir '/' 'aud_vis_loc_lfp'])
cellfun(@(x) copyfile([x '.png'],[cpy_dir '/' 'aud_vis_loc_lfp']),tmp_cpy);

% Find all location files with 'repetition' in lfp
tmp_cpy = strsplit(ls([plt_dir '/' '*' '/' 'initial_location' '/' '*hgp*' '/' '*_rep_*']),'.png');
mkdir([cpy_dir '/' 'aud_vis_loc_hgp'])
cellfun(@(x) copyfile([x '.png'],[cpy_dir '/' 'aud_vis_loc_hgp']),tmp_cpy);