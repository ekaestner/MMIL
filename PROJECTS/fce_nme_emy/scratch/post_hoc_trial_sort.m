clear; clc

dta_loc     = '/home/ekaestner/';
dta_loc_sub = '/Dropbox/McDonald Lab/Erik/Projects/Imaging/fce_nme_emy/';

lst_out = [ dta_loc '/' dta_loc_sub '/' 'Behavior' '/' 'subject_performance' '/' ];

%% Load files
% Subject info
dta_cph = mmil_readtext([ dta_loc '/' dta_loc_sub '/' 'ScanData' '/' 'cipher.tsv'],['\t']);
    dta_cph(cellfun(@isnumeric,dta_cph(:,1)),1) = cellfun(@num2str,dta_cph(cellfun(@isnumeric,dta_cph(:,1)),1),'uni',0);
    
sbj_nme = unique(dta_cph(2:end,1));
    
red_cap_dta = mmil_readtext([ dta_loc '/' dta_loc_sub '/' 'RedcapData' '/' 'RedcapData_pipe.csv'],'|');
    red_cap_dta(cellfun(@isnumeric,red_cap_dta(:,1)),1) = cellfun(@num2str,red_cap_dta(cellfun(@isnumeric,red_cap_dta(:,1)),1),'uni',0);
    
% Encoding trials %%%%%%%%%%%%%
run_one_enc = mmil_readtext([ dta_loc '/' dta_loc_sub '/' 'Behavior' '/' 'run-1_round.tsv' ],['\t']);
run_two_enc = mmil_readtext([ dta_loc '/' dta_loc_sub '/' 'Behavior' '/' 'run-2_round.tsv' ],['\t']);

% Retrieval trials %%%%%%%%%%%


%% Populate behavior
cor_ans = { 'Shawn'  3 ; ...'Shawn'  3 ;
            'Sandra' 2 ; ...'Sandra' 2 ;
            'Eugene' 3 ; ...'Eugene' 3 ;
            'Andrea' 1 ; ...'Connie' 2 ;
            'Daniel' 1 ; ...'Daniel' 1 ;
            'Erica'  2 ; ...'Erica'  2 ;
            'Laura'  3 ; ...'Laura'  3 ;
            'Gerald' 1 ; ...'Joseph' 3 ;
            'Roger'  3 ; ...'Roger'  3 ;
            'Larry'  2 ; ...'Shawn'  1 ;
            'Joseph' 1 ; ...'Daniel' 2 ;
            'Willie' 3 ; ...'Willie' 3 ;
            'Irene'  2 ; ...'Irene'  2 ;
            'Tracy'  3 ; ...'Tracy'  3 ;
            'Connie' 3 ; ...'Joyce'  2 ;
            'Elaine' 2 ; ...'Elaine' 2 ;
            'Monica' 1 ; ...'Joyce'  2 ;
            'Nancy'  2 ; ...'Erica'  1 ;
            'Brenda' 3 ; ...'Carol'  1 ;
            'Dennis' 1 ; ...'Larry'  3 ;
            'Jacob'  2 ; ...'Daniel' 1 ;
            'Harold' 1 ; ...'Harold' 1 ;
            'David'  2 ; ...'David'  2 ;
            'Randy'  1 ; ...'Shawn'  3 ;
            'Scott'  2 ; ...'David'  3 ;
            'Lauren' 1 ; ...'Lauren' 1 ;
            'Wanda'  3 ; ...'Wanda'  3 ;
            'Sheila' 2 ; ...'Elaine' 1 ;
            'Joyce'  3 ; ...'Irene'  2 ;
            'Carol'  1 }; %'Carol'  1 };
        
for iST = 1:size(cor_ans,1)
    
    ind_one = find(strcmpi(run_one_enc(:,5),[cor_ans{iST,1},'.jpg']));
    ind_two = find(strcmpi(run_two_enc(:,5),[cor_ans{iST,1},'.jpg']));
    
    if ~isempty(ind_one)
        cor_ans{iST,3} = 1;
        cor_ans{iST,4} = ind_one;
    elseif ~isempty(ind_two)
        cor_ans{iST,3} = 2;
        cor_ans{iST,4} = ind_two;
    else
        error('missing stimulus')
    end

end

%% Get performance
for iS = 1:numel(sbj_nme)
        
    ind_one = find(strcmpi(red_cap_dta(:,1),sbj_nme{iS}));
    ind_two = find(strcmpi(red_cap_dta(:,1),[sbj_nme{iS} '_2']));
        
    if ~isempty(ind_one)
        
        sbj_run{1} = [ run_one_enc [ {'RecallCorrect'} ; repmat({NaN},size(run_one_enc,1)-1,1) ] [ {'RecognitionCorrect'} ; repmat({NaN},size(run_one_enc,1)-1,1) ] ] ;
        sbj_run{2} = [ run_two_enc [ {'RecallCorrect'} ; repmat({NaN},size(run_two_enc,1)-1,1) ] [ {'RecognitionCorrect'} ; repmat({NaN},size(run_two_enc,1)-1,1) ] ] ;
        
        rcl_ind = string_find(red_cap_dta(1,:),{'free_recall'});
        chc_ind = string_find(red_cap_dta(1,:),{'mult_choice'});
                       
        for iN = 1:numel(rcl_ind)
            if strcmpi(red_cap_dta{ind_one,rcl_ind(iN)},cor_ans{iN,1}) 
                sbj_run{cor_ans{iN,3}}{cor_ans{iN,4},6} = 1;
            else 
               sbj_run{cor_ans{iN,3}}{cor_ans{iN,4},6} = 0;
            end
            if red_cap_dta{ind_one,chc_ind(iN)}==cor_ans{iN,2}
                sbj_run{cor_ans{iN,3}}{cor_ans{iN,4},7} = 1;
            else
                sbj_run{cor_ans{iN,3}}{cor_ans{iN,4},7} = 0;
            end
        end
        
        cell2csv([ lst_out '/' sbj_nme{iS} '_' 'T1' '_' 'run1' '.csv' ],sbj_run{1})
        cell2csv([ lst_out '/' sbj_nme{iS} '_' 'T1' '_' 'run2' '.csv' ],sbj_run{2})
        clear sbj_run
        
    end
                          
    if ~isempty(ind_two)
        
        sbj_run{1} = [ run_one_enc [ {'RecallCorrect'} ; repmat({NaN},size(run_one_enc,1)-1,1) ] [ {'RecognitionCorrect'} ; repmat({NaN},size(run_one_enc,1)-1,1) ] ] ;
        sbj_run{2} = [ run_two_enc [ {'RecallCorrect'} ; repmat({NaN},size(run_two_enc,1)-1,1) ] [ {'RecognitionCorrect'} ; repmat({NaN},size(run_two_enc,1)-1,1) ] ] ;
        
        
        rcl_ind = string_find(red_cap_dta(1,:),{'free_recall'});
        chc_ind = string_find(red_cap_dta(1,:),{'mult_choice'});

        for iN = 1:numel(rcl_ind)
            if strcmpi(red_cap_dta{ind_two,rcl_ind(iN)},cor_ans{iN,1}) 
                sbj_run{cor_ans{iN,3}}{cor_ans{iN,4},6} = 1;
            else
                sbj_run{cor_ans{iN,3}}{cor_ans{iN,4},6} = 0;
            end
            if red_cap_dta{ind_two,chc_ind(iN)}==cor_ans{iN,2}        
                sbj_run{cor_ans{iN,3}}{cor_ans{iN,4},7} = 1;
            else
                sbj_run{cor_ans{iN,3}}{cor_ans{iN,4},7} = 0;
            end
        end
        
        cell2csv([ lst_out '/' sbj_nme{iS} '_' 'T2' '_' 'run1' '.csv' ],sbj_run{1})
        cell2csv([ lst_out '/' sbj_nme{iS} '_' 'T2' '_' 'run2' '.csv' ],sbj_run{2})
        clear sbj_run
                
    end
            
end

%% Save out




