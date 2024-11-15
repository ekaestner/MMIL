

%% Find session directories
for iD = 1:numel(bid_nme)

    ejk_chk_dir([ new_dta_dir '/' bid_nme{iD} ]);

    bid_dir_use = [ bid_dir '/' bid_nme{iD} '/' ];

    sub_dir = dir(bid_dir_use); sub_dir = { sub_dir(:).name }; sub_dir = sub_dir(string_find(sub_dir,'sub-'));

    ses_out = cell(numel(sub_dir),2);
    for iS = 1:numel(sub_dir)
        ses_out{iS,1} = sub_dir{iS};
        ses_dir = dir([ bid_dir_use '/' sub_dir{iS} ]); ses_dir = {ses_dir(:).name};
        if any(string_find(ses_dir,'ses-'))
            ses_out{iS,2} = append(ses_dir(string_find(ses_dir,'ses-')),'|');
            ses_out{iS,2} = [ses_out{iS,2}{:}];
            ses_out{iS,2} = ses_out{iS,2}(1:end-1);
        else
            ses_out{iS,2} = '';
        end
    end

    cell2csv([ new_dta_dir '/' bid_nme{iD} '/' bid_nme{iD} '_' 'sessions' '_' dte_str '.csv'], ses_out);

end

%% Find number of T1's
t1w_num = nan(2,numel(bid_nme));

for iD = 1:numel(bid_nme)

    bid_dir_use = [ bid_dir '/' bid_nme{iD} '/' ];

    ses_dta = mmil_readtext([ new_dta_dir '/' bid_nme{iD} '/' bid_nme{iD} '_' 'sessions' '_' dte_str '.csv']);
    ses_dta = [ ses_dta cell(size(ses_dta,1),6)]; % Adding: # sessions | used session | filenames | # filenames | used filename | errors

    t1w_num(1,iD) = size(ses_dta,1);

    for iS = 1:size(ses_dta,1)

        if rem(iS,100)==0
            fprintf('Dataset: %s | Scan Count: %i from %i\n',bid_nme{iD},iS,size(ses_dta,1))
        end

        if ~isempty(ses_dta{iS,2})
            ses_use_sbj = strsplit(ses_dta{iS,2},'|');
            ses_dta{iS,3} = numel(ses_use_sbj);
            if numel(ses_use_sbj)>1
                ses_use_sbj(string_find(ses_use_sbj,ses_exc.(bid_nme{iD}))) = [];
                if numel(ses_use_sbj)>0
                    ses_use_sbj = ses_use_sbj{1};
                elseif numel(ses_use_sbj)==0
                    ses_use_sbj = '';
                else
                    error("Problem with too many session directories")
                end
            else
                ses_use_sbj = ses_use_sbj{1};
            end
        else
            ses_dta{iS,3} = 0;
            ses_use_sbj = '';
        end

        ses_dta{iS,4} = ses_use_sbj;
        t1w_fle = dir([ bid_dir_use '/' ses_dta{iS,1} '/' ses_use_sbj '/' 'anat' '/' ]); t1w_fle = {t1w_fle(:).name}; t1w_fle(strcmpi(t1w_fle,'.')) = []; t1w_fle(strcmpi(t1w_fle,'..')) = [];
            
            if isempty(t1w_fle); ses_dta{iS,8} = 'Empty anat folder'; continue; end
            if ~isempty(string_find(t1w_fle,'bsub-')); ses_dta{iS,9} = 'CAT12 output??'; end
            
            if isempty(t1w_fle(string_find(t1w_fle,'.nii.gz'))) && ~isempty(t1w_fle(string_find(t1w_fle,'.nii')))
                t1w_fle = t1w_fle(string_find(t1w_fle,'.nii')); ses_dta{iS,8} = '.nii instead of .nii.gz';
            else
                t1w_fle = t1w_fle(string_find(t1w_fle,'.nii.gz'));
            end
            
            ses_dta{iS,5} = append(t1w_fle,'|');
            ses_dta{iS,5} = [ses_dta{iS,5}{:}];
            ses_dta{iS,5} = ses_dta{iS,5}(1:end-1);
            
            ses_dta{iS,6} = numel(t1w_fle);
            
            if isempty(string_find(t1w_fle,'T1w')); ses_dta{iS,8} = 'No identified T1w'; continue; end
            try ses_dta{iS,7} = t1w_fle{string_find(t1w_fle,'T1w')}; catch; error('figure me out')
            end

    end

    cell2csv([ new_dta_dir '/' bid_nme{iD} '/' bid_nme{iD} '_' 'sessions' '_' dte_str '.csv'], ses_dta);
    t1w_num(2,iD) = sum(~cellfun(@isempty,ses_dta(:,7)));

end

t1w_num_out = [ {''} bid_nme {'total'};  {'Folders'} num2cell(t1w_num(1,:)) {sum(t1w_num(1,:))} ; {'T1w_scans'} num2cell(t1w_num(2,:)) {sum(t1w_num(2,:))}];
cell2csv([ new_dta_dir '/' 'T1w' '_' 'count' '_' dte_str '.csv'],t1w_num_out);


