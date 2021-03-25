function out_csv = enigma_relabel(inp_csv)

%% Dx
out_csv(:,1) = inp_csv(:,2);

for iR = 1:size(out_csv,1)
    if out_csv{iR,1}==0
        out_csv{iR,1} = 'Control';
    elseif out_csv{iR,1}==1
        out_csv{iR,1} = 'Epilepsy';
    end
end

%% SDx
out_csv(:,2) = inp_csv(:,3);

for iR = 1:size(out_csv,1)
    if out_csv{iR,2}==0
        out_csv{iR,2} = 'Control';
    elseif out_csv{iR,2}==2
        out_csv{iR,2} = 'GGE';
    elseif out_csv{iR,2}==3
        out_csv{iR,2} = 'L-TLE';
    elseif out_csv{iR,2}==4
        out_csv{iR,2} = 'R-TLE';
    elseif out_csv{iR,2}==5
        out_csv{iR,2} = 'L-TLE';
    elseif out_csv{iR,2}==6
        out_csv{iR,2} = 'R-TLE';
    elseif out_csv{iR,2}==7
        out_csv{iR,2} = 'Extra-Focal';
    elseif out_csv{iR,2}==8
        out_csv{iR,2} = 'Unspecified';
    end
end

%% Sex
out_csv(:,3) = inp_csv(:,6);

for iR = 1:size(out_csv,1)
    if out_csv{iR,3}==1
        out_csv{iR,3} = 'Male';
    elseif out_csv{iR,3}==2
        out_csv{iR,3} = 'Female';
    end
end

%% Handedness
out_csv(:,4) = inp_csv(:,7);
for iR = 1:size(out_csv,1)
    if out_csv{iR,4}==1
        out_csv{iR,4} = 'Right';
    elseif out_csv{iR,4}==2
        out_csv{iR,4} = 'Left';
    end
end


end