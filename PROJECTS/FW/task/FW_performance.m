function FW_performance

clear; clc;
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical';

sbj_nme = mmil_readtext([fcfg.clr_fld '/' 'subjects']);

clear prf_dta lat_dta

for iS = 1:numel(sbj_nme)
    
    log_fle = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' sbj_nme{iS}],'trialf_add');
    
    tar_typ = 0;
    fls_typ = 0;
    rsp_tme = [];
    crr_rsp = 0;
    fls_alr = 0;
    
    log_eve = [];
    
    if ~isempty(log_fle{1})
        for iF = 1:numel(log_fle)
            fwv_trl_hld = mmil_readtext([fcfg.clr_fld '/' '/' 'Logfiles' '/' sbj_nme{iS} '/' log_fle{iF}],'[\t]');
            clear fwv_trl
            if strcmpi(fwv_trl_hld{4,6},'TTime')
                fwv_trl(:,[1 2 4]) = fwv_trl_hld([4:end],[3:4 5]);
            elseif strcmpi(fwv_trl_hld{4,1},'Trial')
                fwv_trl(:,[1 2 4]) = fwv_trl_hld([4:end],[2:3 4]);
            elseif strcmpi(fwv_trl_hld{4,6},'portcode(num)')
                fwv_trl = fwv_trl_hld([4:end],[3:4 6:7]);
            end
            
            rmv_ind = [];
            for iT = 1:size(fwv_trl,1)
                if ~isnumeric(fwv_trl{iT,2}) && ~isempty(string_find(fwv_trl(iT,2),{'mark'}))
                    rmv_ind = [rmv_ind ; iT];
                end
            end
            fwv_trl(rmv_ind,:) = [];
            
            if all(cellfun(@isempty,fwv_trl(:,3))) && strcmpi(fwv_trl_hld{4,6},'TTime')
                for iT = 1:size(fwv_trl,1)
                    if ~isnumeric(fwv_trl{iT,2})
                        if ~cellfun(@isempty,strfind(fwv_trl(iT,2),'T'))
                            fwv_trl{iT,3} = 7;
                        elseif ~cellfun(@isempty,strfind(fwv_trl(iT,2),'O'))
                            fwv_trl{iT,3} = 4;
                        elseif ~cellfun(@isempty,strfind(fwv_trl(iT,2),'npnw'))
                            fwv_trl{iT,3} = 5;
                        elseif ~cellfun(@isempty,strfind(fwv_trl(iT,2),'FF'))
                            fwv_trl{iT,3} = 6;
                        elseif ~cellfun(@isempty,strfind(fwv_trl(iT,2),'nw'))
                            fwv_trl{iT,3} = 3;
                        end
                    elseif ~isempty(fwv_trl{iT,2}) && fwv_trl{iT,2}==64
                        fwv_trl{iT,3} = 64;
                    end
                end
            elseif all(cellfun(@isempty,fwv_trl(:,3))) && strcmpi(fwv_trl_hld{4,1},'Trial')
                for iT = 1:size(fwv_trl,1)
                    if ~isempty(fwv_trl{iT,2}) && ~isnumeric(fwv_trl{iT,2})
                        stm = strsplit(fwv_trl{iT,2},' ');
                        if numel(stm) == 4
                            if ~cellfun(@isempty,strfind(stm(2),'T'))
                                fwv_trl{iT,3} = 7;
                            elseif ~cellfun(@isempty,strfind(stm(2),'O'))
                                fwv_trl{iT,3} = 4;
                            elseif ~cellfun(@isempty,strfind(stm(2),'npnw'))
                                fwv_trl{iT,3} = 5;
                            elseif ~cellfun(@isempty,strfind(stm(2),'w'))
                                fwv_trl{iT,3} = 3;
                            end
                        elseif numel(stm) == 3 && (strcmpi(stm(2),'<<<+>>>') || strcmpi(stm(2),'<<<=>>>'))
                            fwv_trl{iT,3} = 6;
                        end
                    elseif strcmpi(fwv_trl{iT,1},'Response')
                        fwv_trl{iT,2} = 64;
                        fwv_trl{iT,3} = 64;
                    end
                end
            end
                    
            for iT = 2:size(fwv_trl)
                if ~isempty(fwv_trl{iT,2}) && isnumeric(fwv_trl{iT,2}) && fwv_trl{iT,2}==64 && ~isnumeric(fwv_trl{iT-1,2}) && ~strcmpi(fwv_trl(iT-1,2),'break_output') && ~strcmpi(fwv_trl(iT-1,2),'ins trial') && ~strcmpi(fwv_trl(iT-1,2),'break trial') && ~strcmpi(fwv_trl(iT-1,2),'trial_fix')
                    if (~isempty(fwv_trl{iT-2,3}) && fwv_trl{iT-1,3}==7) || (~isempty(fwv_trl{iT-2,3}) && fwv_trl{iT-2,3}==7)
                        crr_rsp = crr_rsp+1;
                        if fwv_trl{iT-1,3}==7
                            rsp_tme = [rsp_tme ; (fwv_trl{iT,4} - fwv_trl{iT-1,4})/10];
                        elseif fwv_trl{iT-2,3}==7
                            rsp_tme = [rsp_tme ; (fwv_trl{iT,4} - fwv_trl{iT-2,4})/10];
                        end
                    else
                        fls_alr = fls_alr+1;
                    end
                elseif fwv_trl{iT,3}==7
                    tar_typ = tar_typ + 1;
                elseif ~isempty(fwv_trl{iT,3})
                    fls_typ = fls_typ + 1;
                end
            end
            
            clear log_eve_hld
            ind = 1;
            for iT = 1:size(fwv_trl,1)
                if isnumeric(fwv_trl{iT,3}) && ~isempty(fwv_trl{iT,3}) && (fwv_trl{iT,3}==3 || fwv_trl{iT,3}==4 || fwv_trl{iT,3}==5 || fwv_trl{iT,3}==6) 
                log_eve_hld(ind) = fwv_trl(iT,3);
                ind = ind + 1;
                end
            end
            log_eve = [log_eve ; cell2mat(log_eve_hld')];
            
        end
        
        rsp_num(iS,1) = crr_rsp;
        rsp_pct(iS,1) = crr_rsp/tar_typ;
        fls_num(iS,1) = fls_alr;
        prf_dta(iS,1) = norminv(crr_rsp/tar_typ-0.0001,0,1) - norminv(fls_alr/fls_typ+0.0001,0,1);
        lat_dta(iS,1) = mean(rsp_tme);
        rsp_dat{iS,1} = rsp_tme;
        
    else
        rsp_num(iS,1) = NaN;
        rsp_pct(iS,1) = NaN;
        fls_num(iS,1) = NaN;
        prf_dta(iS,1) = NaN;
        lat_dta(iS,1) = NaN;
        rsp_dat{iS,1} = [];
        
    end
    
end

save('/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/Logfiles/SubjectPerformance.mat','prf_dta','lat_dta','rsp_dat','rsp_num','rsp_pct','fls_num')

% Plot
figure('Visible','off')

subplot(1,3,1)
hold on
for iS = 1:size(prf_dta,1) 
    scatter((1.1-0.9).*rand(1)+0.9,prf_dta(iS,1),80,'filled','MarkerEdgeColor','k','MarkerFaceColor','k');
end
xlim([0.8 1.2]); set(gca,'XColor',[1 1 1]);
ylim([0 5.5]); ylabel('D''');

subplot(1,3,2)
hold on
for iS = 1:size(lat_dta,1) 
    xxx_val = (1.1-0.9).*rand(1)+0.9;
    scatter(xxx_val,lat_dta(iS,1),80,'filled','MarkerEdgeColor','k','MarkerFaceColor','k');
    if ~isnan(lat_dta(iS,1))
        line([xxx_val xxx_val],[lat_dta(iS,1)-std(rsp_dat{iS,1})/sqrt(numel(rsp_dat{iS,1})) lat_dta(iS,1)+std(rsp_dat{iS,1})/sqrt(numel(rsp_dat{iS,1}))],'Color','k');
        line([xxx_val-0.0015 xxx_val+0.0015],[lat_dta(iS,1)+std(rsp_dat{iS,1})/sqrt(numel(rsp_dat{iS,1})) lat_dta(iS,1)+std(rsp_dat{iS,1})/sqrt(numel(rsp_dat{iS,1}))],'Color','k');
        line([xxx_val-0.0015 xxx_val+0.0015],[lat_dta(iS,1)-std(rsp_dat{iS,1})/sqrt(numel(rsp_dat{iS,1})) lat_dta(iS,1)-std(rsp_dat{iS,1})/sqrt(numel(rsp_dat{iS,1}))],'Color','k');
    end
end
xlim([0.8 1.2]); set(gca,'XColor',[1 1 1]);
ylim([0 1200]); ylabel('Response Time (ms)');

tightfig()

pcfg = [];
pcfg.fle_nme = ['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/figure1/FWbehavior'];
pcfg.prn_typ = 'eps';
pcfg.jpg     = 0;
pcfg.cst_col = 0;
mmil_print_plot(pcfg)

pcfg = [];
pcfg.fle_nme = ['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/figure1/FWbehavior'];
pcfg.prn_typ = 'png';
pcfg.jpg     = 0;
pcfg.cst_col = 0;
mmil_print_plot(pcfg)

close all

% Save output
for iPD = 1:numel(prf_dta)
prf_dta_cll{iPD} = sprintf('%g',roundsd(prf_dta(iPD),3));
end

beh_dta = [ sbj_nme ...
      prf_dta_cll'...
      num2cell(round(lat_dta)) ...
      cellfun(@(x) [num2str(x) '%'],num2cell(round(rsp_pct*100)),'uni',0) ...
      num2cell(rsp_num) num2cell(fls_num)];

  for iR = 1:size(beh_dta,1)
      for iC = 1:size(beh_dta,2)
          if ischar(beh_dta{iR,iC})
              if strcmpi(beh_dta{iR,iC},'NaN%') || strcmpi(beh_dta{iR,iC},'NaN') || strcmpi(beh_dta{iR,iC},'0%')
                   beh_dta{iR,iC} = '-';
              end
          elseif isnumeric(beh_dta{iR,iC})
              if isnan(beh_dta{iR,iC}) || ((iC == 4 || iC == 5) && beh_dta{iR,iC}<1 )
                  beh_dta{iR,iC} = '-';
              end
          end
      end
  end
      
cell2csv([fcfg.clr_fld '/' 'manuscript' '/' 'supplmental_tables' '/' 'behavior.csv'], ...
         beh_dta )



end

