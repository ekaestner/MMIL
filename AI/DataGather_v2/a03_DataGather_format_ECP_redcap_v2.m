
%% Missing
red_cap_col = [ red_cap_col {'duration_add'} {'localization_add'} {'focal_add'} {'img_nme_add'} ];
red_cap_dta = [ red_cap_dta cell(size(red_cap_dta,1),4) ];

%% Calculations
hnd_col = strcmpi(red_cap_col,'ehq');

dur_col = strcmpi(red_cap_col,'duration_add');
age_col = strcmpi(red_cap_col,'age');
ons_col = strcmpi(red_cap_col,'ao1');

dgn_col = strcmpi(red_cap_col,'dx');
loc_col = strcmpi(red_cap_col,'localization_add');
foc_col = strcmpi(red_cap_col,'focal_add');

nme_col = strcmpi(red_cap_col,'img_nme_add');

for iS = 1:numel(red_cap_sbj)

    %% Calculate Handedness
    if ~isempty(red_cap_dta{iS,hnd_col})
        if red_cap_dta{iS,hnd_col}>40
            red_cap_dta{iS,hnd_col} = 'Right';
        elseif red_cap_dta{iS,hnd_col}<-40
            red_cap_dta{iS,hnd_col} = 'Left';
        else
            red_cap_dta{iS,hnd_col} = 'Ambi';
        end
    else
        red_cap_dta{iS,hnd_col} = '';
    end

    %% Calculate Duration
    if ~isempty(red_cap_dta{iS,age_col}) && ~isempty(red_cap_dta{iS,ons_col})
        red_cap_dta{iS,dur_col} = red_cap_dta{iS,age_col} - red_cap_dta{iS,ons_col};
    end

    %% Add Focality & Localization
    if ~isempty(red_cap_dta{iS,dgn_col}) && red_cap_dta{iS,dgn_col}==1
        red_cap_dta{iS,foc_col} = 'focal';
        red_cap_dta{iS,loc_col} = 'temporal';
    else
        red_cap_dta{iS,loc_col} = '';
        red_cap_dta{iS,foc_col} = '';
    end

    %% Imaging Name
    red_cap_dta{iS,nme_col} = [ 'ECP' red_cap_sbj{iS,1} ];

end








