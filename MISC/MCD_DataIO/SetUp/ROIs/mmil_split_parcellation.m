% File              : mmil_split_parcellation.m
% Author            : Erik Kaestner <erik.kaestner@gmail.com>
% Date              : 08.11.2018
% Last Modified Date: 15.11.2018
% Last Modified By  : Akshara Balachandra <abalacha@ucsd.edu>

function mmil_split_parcellation(cfg)

if ~isfield( cfg , 'ovr_wrt' ); cfg.ovr_wrt = 0; end

hms = { 'lh' 'rh' };

fprintf( [ cfg.sbj_nme ' : ' 'Starting ' ] )

if ~( exist([ cfg.prj_dir '/' 'DATA' '/' cfg.sbj_nme '/' 'ROIs' '/' 'lh.aparc.split.annot' ], 'file') == 2 ) || cfg.ovr_wrt == 1
    
    if cfg.ovr_wrt == 1 && exist([ cfg.prj_dir '/' 'DATA' '/' cfg.sbj_nme '/' 'ROIs' '/' 'lh.aparc.split.annot' ], 'file') 
        delete( [ cfg.prj_dir '/' 'DATA' '/' cfg.sbj_nme '/' 'ROIs' '/' 'lh.aparc.split.annot' ] )
        delete( [ cfg.prj_dir '/' 'DATA' '/' cfg.sbj_nme '/' 'ROIs' '/' 'rh.aparc.split.annot' ] )
    end

    if ~strcmpi( cfg.sbj_nme , 'fsaverage' )
        sbj_dir_lst = dir(sprintf('%s/FSURF_*',cfg.fsr_dir));
        sbj_dir_lst = regexp({sbj_dir_lst.name},['FSURF_' cfg.fsr_nme '.+_1$'],'match'); sbj_dir_lst = [sbj_dir_lst{:}];
    else
        sbj_dir_lst = dir(sprintf('%s/*',cfg.fsr_dir));
        sbj_dir_lst = regexp({sbj_dir_lst.name},['' cfg.fsr_nme],'match'); sbj_dir_lst = [sbj_dir_lst{:}];
    end
    
    if ~isempty(sbj_dir_lst)
        
        for iH = 1:length(hms)
            
            %% Make Split
            out_fle = [ cfg.prj_dir '/' 'DATA' '/' cfg.sbj_nme '/' 'ROIs' '/' hms{iH} '.' 'aparc' '.' 'split' '.' 'annot'];
            spl_fle = '/home/ekaestne/PROJECTS/DATA/mmil_split_cfg.txt';
            
            %% csh
            %% cmd = sprintf('setenv SUBJECTS_DIR %s \n',out_pth);
            %% bash/zsh/sh
            
            cmd = sprintf( 'setenv SUBJECTS_DIR %s \n' , cfg.fsr_dir );
            cmd = [ cmd 'mris_divide_parcellation' ' ' sbj_dir_lst{1} ' ' hms{iH} ' ' 'aparc.annot' ' '...
                spl_fle ' ' out_fle ' ' '--sbjdir' ' ' cfg.fsr_dir ];
            [status,msg] = unix(cmd);
            if status; error(msg); end
            
            %% change colortable.struct_names
            [vertices, label, colortable] = read_annotation(out_fle);
            
            colortable.struct_names = regexprep(colortable.struct_names,'superiortemporal_div1','caudal-STG');
            colortable.struct_names = regexprep(colortable.struct_names,'superiortemporal_div2','middle-STG');
            colortable.struct_names = regexprep(colortable.struct_names,'superiortemporal_div3','rostral-STG');
            
            colortable.struct_names = regexprep(colortable.struct_names,'middletemporal_div1','caudal-MTG');
            colortable.struct_names = regexprep(colortable.struct_names,'middletemporal_div2','middle-MTG');
            colortable.struct_names = regexprep(colortable.struct_names,'middletemporal_div3','rostral-MTG');
            
            colortable.struct_names = regexprep(colortable.struct_names,'inferiortemporal_div1','caudal-ITG');
            colortable.struct_names = regexprep(colortable.struct_names,'inferiortemporal_div2','middle-ITG');
            colortable.struct_names = regexprep(colortable.struct_names,'inferiortemporal_div3','rostral-ITG');
            
            colortable.struct_names = regexprep(colortable.struct_names,'fusiform_div1','caudal-fusiform');
            colortable.struct_names = regexprep(colortable.struct_names,'fusiform_div2','middle-fusiform');
            colortable.struct_names = regexprep(colortable.struct_names,'fusiform_div3','rostral-fusiform');
            
            colortable.struct_names = regexprep(colortable.struct_names,'precentral_div1','superior-precentral');
            colortable.struct_names = regexprep(colortable.struct_names,'precentral_div2','middle-precentral');
            colortable.struct_names = regexprep(colortable.struct_names,'precentral_div3','inferior-precentral');
            
            colortable.struct_names = regexprep(colortable.struct_names,'postcentral_div1','superior-postcentral');
            colortable.struct_names = regexprep(colortable.struct_names,'postcentral_div2','middle-postcentral');
            colortable.struct_names = regexprep(colortable.struct_names,'postcentral_div3','inferior-postcentral');
            
            colortable.struct_names = regexprep(colortable.struct_names,'rostralmiddlefrontal_div1','middle-middlefrontal');
            colortable.struct_names = regexprep(colortable.struct_names,'rostralmiddlefrontal_div2','rostral-middlefrontal');
            
%             colortable.struct_names = regexprep(colortable.struct_names,'parahippocampal_div1','caudal-parahippocampal');
%             colortable.struct_names = regexprep(colortable.struct_names,'parahippocampal_div2','middle-parahippocampal');
            
            colortable.struct_names = regexprep(colortable.struct_names,'superiorfrontal_div1','caudal-superiorfrontal');
            colortable.struct_names = regexprep(colortable.struct_names,'superiorfrontal_div2','middle-superiorfrontal');
            colortable.struct_names = regexprep(colortable.struct_names,'superiorfrontal_div3','rostral-superiorfrontal');
            
            %% change colortable.table
            idx1 = find(strcmp(colortable.struct_names,'middle-STG'));
            mSTG_id = colortable.table(idx1,5);
            colortable.table(idx1,1:3) = [24,200,24];
            colortable.table(idx1,5) = colortable.table(idx1,1:4)*[1;2^8;2^16;2^24];
            
            idx2 = find(strcmp(colortable.struct_names,'middle-MTG'));
            mMTG_id = colortable.table(idx2,5);
            colortable.table(idx2,1:3) = [110,160,220];
            colortable.table(idx2,5) = colortable.table(idx2,1:4)*[1;2^8;2^16;2^24];
            
            idx3 = find(strcmp(colortable.struct_names,'middle-ITG'));
            mITG_id = colortable.table(idx3,5);
            colortable.table(idx3,1:3) = [75,0,110];
            colortable.table(idx3,5) = colortable.table(idx3,1:4)*[1;2^8;2^16;2^24];
            
            idx4 = find(strcmp(colortable.struct_names,'middle-fusiform'));
            mfusiform_id = colortable.table(idx4,5);
            colortable.table(idx4,1:3) = [21,68,6];
            colortable.table(idx4,5) = colortable.table(idx4,1:4)*[1;2^8;2^16;2^24];
            
            idx5 = find(strcmp(colortable.struct_names,'middle-precentral'));
            mprecentral_id = colortable.table(idx5,5);
            colortable.table(idx5,1:3) = [255,255,20];
            colortable.table(idx5,5) = colortable.table(idx5,1:4)*[1;2^8;2^16;2^24];
            
            idx6 = find(strcmp(colortable.struct_names,'middle-postcentral'));
            mpostcentral_id = colortable.table(idx6,5);
            colortable.table(idx6,1:3) = [173,10,253];
            colortable.table(idx6,5) = colortable.table(idx6,1:4)*[1;2^8;2^16;2^24];
            
            idx7 = find(strcmp(colortable.struct_names,'middle-middlefrontal'));
            mmiddlefrontal_id = colortable.table(idx7,5);
            colortable.table(idx7,1:3) = [10,136,138];
            colortable.table(idx7,5) = colortable.table(idx7,1:4)*[1;2^8;2^16;2^24];
            
            idx8 = find(strcmp(colortable.struct_names,'middle-superiorfrontal'));
            msuperiorfrontal_id = colortable.table(idx8,5);
            colortable.table(idx8,1:3) = [190, 3,253];
            colortable.table(idx8,5) = colortable.table(idx8,1:4)*[1;2^8;2^16;2^24];
            
%             idx9 = find(strcmp(colortable.struct_names,'caudal-parahippocampal'));
%             mmiddleparahip_id = colortable.table(idx9,5);
%             colortable.table(idx9,1:3) = [201,10,130];
%             colortable.table(idx9,5) = colortable.table(idx9,1:4)*[1;2^8;2^16;2^24];
%             
%             idx10 = find(strcmp(colortable.struct_names,'middle-parahippocampal'));
%             caudparahip_id = colortable.table(idx10,5);
%             colortable.table(idx10,1:3) = [10,144,253];
%             colortable.table(idx10,5) = colortable.table(idx10,1:4)*[1;2^8;2^16;2^24];
            
            %% replace ID in label
            label(label==mSTG_id) = colortable.table(idx1,5);
            label(label==mMTG_id) = colortable.table(idx2,5);
            label(label==mITG_id) = colortable.table(idx3,5);
            label(label==mfusiform_id)      = colortable.table(idx4,5);
            label(label==mprecentral_id)    = colortable.table(idx5,5);
            label(label==mpostcentral_id)   = colortable.table(idx6,5);
            label(label==mmiddlefrontal_id) = colortable.table(idx7,5);
            label(label==msuperiorfrontal_id) = colortable.table(idx8,5);
%             label(label==mmiddleparahip_id) = colortable.table(idx9,5);
%             label(label==caudparahip_id) = colortable.table(idx10,5);
            
            %% replace ID in vertices - Not necessary
            %         vertices(vertices==mSTG_id) = colortable.table(idx1,5);
            %         vertices(vertices==mMTG_id) = colortable.table(idx2,5);
            %         vertices(vertices==mITG_id) = colortable.table(idx3,5);
            %         vertices(vertices==mfusiform_id)    = colortable.table(idx4,5);
            %         vertices(vertices==mprecentral_id)  = colortable.table(idx5,5);
            %         vertices(vertices==mpostcentral_id) = colortable.table(idx6,5);
            %         vertices(vertices==mmiddlefrontal_id) = colortable.table(idx7,5);
            %         vertices(vertices==msuperiorfrontal_id) = colortable.table(idx8,5);
            %         vertices(vertices==mmiddleparahip_id) = colortable.table(idx9,5);
            %         vertices(vertices==caudparahip_id) = colortable.table(idx10,5);
            
            %% save new annotation
            colortable.numEntries   = 59;
            colortable.struct_names = colortable.struct_names(1:59);
            colortable.table        = colortable.table(1:59,:);
            
            write_annotation(out_fle,vertices,label,colortable);
            
            fprintf(' - finished! \n')
            
        end
    else
        fprintf(' - %%%% MISSING %%%% \n')
    end
else
    
    fprintf(' - Did not perform \n')
    
end

end

