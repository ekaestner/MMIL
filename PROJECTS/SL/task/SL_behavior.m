function SL_behavior

clear; clc;
fcfg.clr_fld = '/home/ekaestne/PROJECTS/OUTPUT/SL';

sbj_nme = mmil_readtext([fcfg.clr_fld '/' 'subjects']);

clear prf_dta lat_dta rsp_dat

prf_dta = cell(1,3);
lat_dta = cell(1,3);
rsp_dat = cell(1,3);

for iS = 1:numel(sbj_nme)
   
    log_fle = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' sbj_nme{iS}],'trialf_add');
    
    tar_typ = [0 0 0 0];
    fls_typ = [0 0 0 0];
    rsp_tme = {[] [] [] []};
    crr_rsp = [0 0 0 0];
    fls_alr = [0 0 0 0];
    
    log_eve = [];
    
    if ~isempty(log_fle{1})
        for iF = 1:numel(log_fle)
                       
            fwv_trl_hld = mmil_readtext([fcfg.clr_fld '/' '/' 'Logfiles' '/' sbj_nme{iS} '/' mmil_spec_char(sbj_nme{iS},{'_'},{'-'}) '_Vis1st_' log_fle{iF} '.log'],'[\t]');
            clear fwv_trl
            if strcmpi(fwv_trl_hld{4,6},'TTime')
                fwv_trl(:,[1 2 4]) = fwv_trl_hld([4:end],[3:4 5]);
            elseif strcmpi(fwv_trl_hld{4,1},'Trial')
                fwv_trl(:,[1 2 4]) = fwv_trl_hld([4:end],[2:3 4]);
            elseif strcmpi(fwv_trl_hld{4,6},'portcode(num)')
                fwv_trl = fwv_trl_hld([4:end],[3:4 6:7]);
            end
            
            fst_one = find(cellfun(@isstr,fwv_trl(:,2)));            
            fwv_trl = fwv_trl(fst_one(2)-2:end,:);
            fst_one = find(cellfun(@(x) isnumeric(x) && (x==1 || x==2 || x==3 || x==4),fwv_trl(:,2))) ;
            
            for iT = 1:numel(fst_one)
                                
                % Match/Mismatch
                if fwv_trl{fst_one(iT),2}==1
                    tar_typ(1) = tar_typ(1)+1;
                    fls_typ(2) = fls_typ(2)+1;
                    if strcmpi(fwv_trl{fst_one(iT)+2,2},'correct')
                        crr_rsp(1) = crr_rsp(1)+1;
                        rsp_tme{1} = [rsp_tme{1} fwv_trl{fst_one(iT)+1,4}/10-fwv_trl{fst_one(iT),4}/10-450];
                    else
                        fls_alr(1) = fls_alr(1)+1; 
                    end
                    
                    %
                elseif fwv_trl{fst_one(iT),2}==2
                    tar_typ(2) = tar_typ(2)+1;
                    fls_typ(1) = fls_typ(1)+1;
                    if strcmpi(fwv_trl{fst_one(iT)+2,2},'correct')
                        crr_rsp(2) = crr_rsp(2)+1;
                        rsp_tme{2} = [rsp_tme{2} fwv_trl{fst_one(iT)+1,4}/10-fwv_trl{fst_one(iT),4}/10-450];
                    else
                        fls_alr(2) = fls_alr(2)+1; 
                    end
                    
                    % Visual Control
                elseif fwv_trl{fst_one(iT),2}==3
                    tar_typ(3) = tar_typ(3)+1;
                    if strcmpi(fwv_trl{fst_one(iT)+2,2},'correct')
                        crr_rsp(3) = crr_rsp(3)+1;
                        rsp_tme{3} = [rsp_tme{3} fwv_trl{fst_one(iT)+1,4}/10-fwv_trl{fst_one(iT),4}/10-450];
                    else
                        fls_alr(3) = fls_alr(3)+1; 
                    end
                    
                    % Auditory Control
                elseif fwv_trl{fst_one(iT),2}==4
                    tar_typ(4) = tar_typ(4)+1;
                    if strcmpi(fwv_trl{fst_one(iT)+2,2},'correct')
                        crr_rsp(4) = crr_rsp(4)+1;
                        rsp_tme{4} = [rsp_tme{4} fwv_trl{fst_one(iT)+1,4}/10-fwv_trl{fst_one(iT),4}/10-450];
                    else
                        fls_alr(4) = fls_alr(4)+1; 
                    end
                    
                end
                
            end
            
        end
        
        prf_dta{1}(iS,1) = norminv(crr_rsp(1)/tar_typ(1)-0.0001,0,1) - norminv(fls_alr(1)/fls_typ(1)+0.0001,0,1);
        prp_dta{1}(iS,1) = crr_rsp(1)/tar_typ(1);
        lat_dta{1}(iS,1) = mean(rsp_tme{1});
        rsp_dat{1}{iS,1} = rsp_tme{1};
        
        prf_dta{2}(iS,1) = norminv(crr_rsp(2)/tar_typ(2)-0.0001,0,1) - norminv(fls_alr(2)/fls_typ(2)+0.0001,0,1);
        prp_dta{2}(iS,1) = crr_rsp(2)/tar_typ(2);
        lat_dta{2}(iS,1) = mean(rsp_tme{2});
        rsp_dat{2}{iS,1} = rsp_tme{2};
        
        prf_dta{3}(iS,1) = norminv((crr_rsp(1)+crr_rsp(2))/(tar_typ(1)+tar_typ(2))-0.0001,0,1) - norminv((fls_alr(1)+fls_alr(2))/(fls_typ(1)+fls_typ(2))+0.0001,0,1);
        prp_dta{3}(iS,1) = (crr_rsp(1)+crr_rsp(2))/(tar_typ(1)+tar_typ(2));
        lat_dta{3}(iS,1) = mean([rsp_tme{1} rsp_tme{2}]);
        rsp_dat{3}{iS,1} = [rsp_tme{1} rsp_tme{2}];
        
        prf_dta{4}(iS,1) = norminv(crr_rsp(3)/tar_typ(3)-0.0001,0,1) - norminv(fls_alr(3)/fls_typ(1)+0.0001,0,1);
        prp_dta{4}(iS,1) = crr_rsp(3)/tar_typ(3);
        lat_dta{4}(iS,1) = mean(rsp_tme{3});
        rsp_dat{4}{iS,1} = rsp_tme{3};
        
        prf_dta{5}(iS,1) = norminv(crr_rsp(4)/tar_typ(4)-0.0001,0,1) - norminv(fls_alr(4)/fls_typ(1)+0.0001,0,1);
        prp_dta{5}(iS,1) = crr_rsp(4)/tar_typ(4);
        lat_dta{5}(iS,1) = mean(rsp_tme{4});
        rsp_dat{5}{iS,1} = rsp_tme{4};
        
    else
        error('')
    end
       
end

%% Plot
% d'
figure('Visible','off')
subplot(1,3,1)
hold on
for iS = 1:size(prf_dta{1},1) 
    scatter((1.1-0.9).*rand(1)+0.9,prf_dta{3}(iS,1),80,'filled','MarkerEdgeColor','k','MarkerFaceColor','k');
    scatter((1.6-1.4).*rand(1)+1.4,prf_dta{4}(iS,1),80,'filled','MarkerEdgeColor','k','MarkerFaceColor','m');
    scatter((2.1-1.9).*rand(1)+1.9,prf_dta{5}(iS,1),80,'filled','MarkerEdgeColor','k','MarkerFaceColor','c');
end
xlim([0.8 2.2]); set(gca,'XColor',[1 1 1]);
ylim([0 5.5]); ylabel('D''');

% Proportion Correct
figure('Visible','off')
subplot(1,3,1)
hold on
for iS = 1:size(prf_dta{1},1) 
    scatter((1.1-0.9).*rand(1)+0.9,prp_dta{1}(iS,1)*100,80,'filled','MarkerEdgeColor','k','MarkerFaceColor','g');
    scatter((1.6-1.4).*rand(1)+1.4,prp_dta{2}(iS,1)*100,80,'filled','MarkerEdgeColor','k','MarkerFaceColor','y');
end
xlim([0.8 1.7]); set(gca,'XColor',[1 1 1]);
ylim([0 100]); ylabel('Correct %');

pcfg = [];
pcfg.fle_nme = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/manuscript/figure1/SLbehavior_proportion'];
pcfg.prn_typ = 'eps';
pcfg.jpg     = 0;
pcfg.cst_col = 0;
mmil_print_plot(pcfg);

% Reaction Time - proportion
figure('Visible','off')
subplot(1,3,1)
hold on
for iS = 1:size(lat_dta{3},1) 
    xx1_val = (1.05-0.95).*rand(1)+0.9;
    xx2_val = (1.55-1.45).*rand(1)+1.4;
        
    if ~isnan(lat_dta{1}(iS,1))
        line([xx1_val xx1_val],[lat_dta{1}(iS,1)-std(rsp_dat{1}{iS,1})/sqrt(numel(rsp_dat{1}{iS,1})) lat_dta{1}(iS,1)+std(rsp_dat{1}{iS,1})/sqrt(numel(rsp_dat{1}{iS,1}))],'Color','k');
        line([xx1_val-0.0015 xx1_val+0.0015],[lat_dta{1}(iS,1)+std(rsp_dat{1}{iS,1})/sqrt(numel(rsp_dat{1}{iS,1})) lat_dta{1}(iS,1)+std(rsp_dat{1}{iS,1})/sqrt(numel(rsp_dat{1}{iS,1}))],'Color','k');
        line([xx1_val-0.0015 xx1_val+0.0015],[lat_dta{1}(iS,1)-std(rsp_dat{1}{iS,1})/sqrt(numel(rsp_dat{1}{iS,1})) lat_dta{1}(iS,1)-std(rsp_dat{1}{iS,1})/sqrt(numel(rsp_dat{1}{iS,1}))],'Color','k');
    end
    scatter(xx1_val,lat_dta{1}(iS,1),80,'filled','MarkerEdgeColor','k','MarkerFaceColor','g');
    
    if ~isnan(lat_dta{2}(iS,1))
        line([xx2_val xx2_val],[lat_dta{2}(iS,1)-std(rsp_dat{2}{iS,1})/sqrt(numel(rsp_dat{2}{iS,1})) lat_dta{2}(iS,1)+std(rsp_dat{2}{iS,1})/sqrt(numel(rsp_dat{2}{iS,1}))],'Color','k');
        line([xx2_val-0.0015 xx2_val+0.0015],[lat_dta{2}(iS,1)+std(rsp_dat{2}{iS,1})/sqrt(numel(rsp_dat{2}{iS,1})) lat_dta{2}(iS,1)+std(rsp_dat{2}{iS,1})/sqrt(numel(rsp_dat{2}{iS,1}))],'Color','k');
        line([xx2_val-0.0015 xx2_val+0.0015],[lat_dta{2}(iS,1)-std(rsp_dat{2}{iS,1})/sqrt(numel(rsp_dat{2}{iS,1})) lat_dta{2}(iS,1)-std(rsp_dat{2}{iS,1})/sqrt(numel(rsp_dat{2}{iS,1}))],'Color','k');
    end
    scatter(xx2_val,lat_dta{2}(iS,1),80,'filled','MarkerEdgeColor','k','MarkerFaceColor','y');
    
end
xlim([0.8 1.7]); set(gca,'XColor',[1 1 1]);
ylim([0 1000]); ylabel('Response Time (ms)');

pcfg = [];
pcfg.fle_nme = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/manuscript/figure1/SLbehavior_responsetime'];
pcfg.prn_typ = 'eps';
pcfg.jpg     = 0;
pcfg.cst_col = 0;
mmil_print_plot(pcfg);

% Reaction Time - all
figure('Visible','off')
subplot(1,3,1)
hold on
for iS = 1:size(lat_dta{3},1) 
    xx1_val = (1.1-0.9).*rand(1)+0.9;
    xx2_val = (1.6-1.4).*rand(1)+1.4;
    xx3_val = (2.1-1.9).*rand(1)+1.9;
    
    scatter(xx1_val,lat_dta{3}(iS,1),80,'filled','MarkerEdgeColor','k','MarkerFaceColor','k');
    if ~isnan(lat_dta{3}(iS,1))
        line([xx1_val xx1_val],[lat_dta{3}(iS,1)-std(rsp_dat{3}{iS,1})/sqrt(numel(rsp_dat{3}{iS,1})) lat_dta{3}(iS,1)+std(rsp_dat{3}{iS,1})/sqrt(numel(rsp_dat{3}{iS,1}))],'Color','k');
        line([xx1_val-0.0015 xx1_val+0.0015],[lat_dta{3}(iS,1)+std(rsp_dat{3}{iS,1})/sqrt(numel(rsp_dat{3}{iS,1})) lat_dta{3}(iS,1)+std(rsp_dat{3}{iS,1})/sqrt(numel(rsp_dat{3}{iS,1}))],'Color','k');
        line([xx1_val-0.0015 xx1_val+0.0015],[lat_dta{3}(iS,1)-std(rsp_dat{3}{iS,1})/sqrt(numel(rsp_dat{3}{iS,1})) lat_dta{3}(iS,1)-std(rsp_dat{3}{iS,1})/sqrt(numel(rsp_dat{3}{iS,1}))],'Color','k');
    end
    
    scatter(xx2_val,lat_dta{4}(iS,1),80,'filled','MarkerEdgeColor','m','MarkerFaceColor','m');
    if ~isnan(lat_dta{4}(iS,1))
        line([xx2_val xx2_val],[lat_dta{4}(iS,1)-std(rsp_dat{4}{iS,1})/sqrt(numel(rsp_dat{4}{iS,1})) lat_dta{4}(iS,1)+std(rsp_dat{4}{iS,1})/sqrt(numel(rsp_dat{4}{iS,1}))],'Color','m');
        line([xx2_val-0.0015 xx2_val+0.0015],[lat_dta{4}(iS,1)+std(rsp_dat{4}{iS,1})/sqrt(numel(rsp_dat{4}{iS,1})) lat_dta{4}(iS,1)+std(rsp_dat{4}{iS,1})/sqrt(numel(rsp_dat{4}{iS,1}))],'Color','m');
        line([xx2_val-0.0015 xx2_val+0.0015],[lat_dta{4}(iS,1)-std(rsp_dat{4}{iS,1})/sqrt(numel(rsp_dat{4}{iS,1})) lat_dta{4}(iS,1)-std(rsp_dat{4}{iS,1})/sqrt(numel(rsp_dat{4}{iS,1}))],'Color','m');
    end
    
    scatter(xx3_val,lat_dta{5}(iS,1),80,'filled','MarkerEdgeColor','c','MarkerFaceColor','c');
    if ~isnan(lat_dta{5}(iS,1))
        line([xx3_val xx3_val],[lat_dta{5}(iS,1)-std(rsp_dat{5}{iS,1})/sqrt(numel(rsp_dat{5}{iS,1})) lat_dta{5}(iS,1)+std(rsp_dat{5}{iS,1})/sqrt(numel(rsp_dat{5}{iS,1}))],'Color','c');
        line([xx3_val-0.0015 xx3_val+0.0015],[lat_dta{5}(iS,1)+std(rsp_dat{5}{iS,1})/sqrt(numel(rsp_dat{5}{iS,1})) lat_dta{5}(iS,1)+std(rsp_dat{5}{iS,1})/sqrt(numel(rsp_dat{5}{iS,1}))],'Color','c');
        line([xx3_val-0.0015 xx3_val+0.0015],[lat_dta{5}(iS,1)-std(rsp_dat{5}{iS,1})/sqrt(numel(rsp_dat{5}{iS,1})) lat_dta{5}(iS,1)-std(rsp_dat{5}{iS,1})/sqrt(numel(rsp_dat{5}{iS,1}))],'Color','c');
    end
    
end
xlim([0.8 2.2]); set(gca,'XColor',[1 1 1]);
ylim([0 1000]); ylabel('Response Time (ms)');

tightfig()

pcfg = [];
pcfg.fle_nme = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/manuscript/figure1/SLbehavior'];
pcfg.prn_typ = 'eps';
pcfg.jpg     = 0;
pcfg.cst_col = 0;
mmil_print_plot(pcfg);

pcfg = [];
pcfg.fle_nme = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/manuscript/figure1/SLbehavior'];
pcfg.prn_typ = 'png';
pcfg.jpg     = 0;
pcfg.cst_col = 0;
mmil_print_plot(pcfg);

close all

%% Tables
beh_dta = [ sbj_nme ...
            cellfun(@(x) [num2str(x) '%'],num2cell(round(prp_dta{1}*100)),'uni',0)...
            num2cell(round(lat_dta{1})) ...
            cellfun(@(x) [num2str(x) '%'],num2cell(round(prp_dta{2}*100)),'uni',0)...
            num2cell(round(lat_dta{2})) ...
            cellfun(@(x) [num2str(x) '%'],num2cell(round(prp_dta{4}*100)),'uni',0)...
            num2cell(round(lat_dta{4})) ...
            cellfun(@(x) [num2str(x) '%'],num2cell(round(prp_dta{5}*100)),'uni',0)...
            num2cell(round(lat_dta{5}))];
 
cell2csv([fcfg.clr_fld '/' 'manuscript' '/' 'supplmental_tables' '/' 'behavior.csv'], ...
         beh_dta )

beh_dta = [ sbj_nme ...
            num2cell(round(prp_dta{1}*100)) ...
            num2cell(round(lat_dta{1})) ...
            num2cell(round(prp_dta{2}*100))...
            num2cell(round(lat_dta{2})) ...
           num2cell(round(prp_dta{4}*100))...
            num2cell(round(lat_dta{4})) ...
            num2cell(round(prp_dta{5}*100))...
            num2cell(round(lat_dta{5}))];
 
cell2csv([fcfg.clr_fld '/' 'manuscript' '/' 'supplmental_tables' '/' 'behavior_num.csv'], ...
         beh_dta )
     
end


