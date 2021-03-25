function [ ] = write_dicoms(img,sl, slice_demux_coord, dataFileHeader, t_array, num_t2, dcmfile, maxval, new_seuid)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    %%%setting DICOM dictionary to GE dicom dictionary
    dicomdict('set','gems-dicom-dict.txt');
    metadata=dicominfo(dcmfile);
    mux_factor=dataFileHeader.rdb_hdr.user6;
    diff_dir = metadata.UserData24;
    
    nslices=(dataFileHeader.rdb_hdr.nslices/dataFileHeader.rdb_hdr.reps)*mux_factor;
    num_rows = size(img,1);
    num_cols = size(img,2);
    
    pixel_spacingx = metadata.DisplayFieldOfView/num_rows;
    pixel_spacingy = metadata.DisplayFieldOfView/num_cols;
    pixel_spacing=[pixel_spacingx;pixel_spacingy];
    pixel_bw = (metadata.PixelBandwidth*num_rows)/256;

    total_slices = nslices;
    total_images = total_slices*(length(t_array)); %MJM
    if strcmp(dataFileHeader.image.psd_iname, 'EPI2')
        
        diffu_dir = metadata.UserData24;
        if dataFileHeader.rdb_hdr.user1 > 0
            max_bval = dataFileHeader.rdb_hdr.user1;
        else
            if diffu_dir < 25
                max_bval = 1000;
            else
                max_bval = 2800;
            end
        end
    end

    outname=sprintf('E%dS%dI',dataFileHeader.exam.ex_no, dataFileHeader.series.se_no*100);

    % Scale pixel intensity
    img=int16((img/maxval)*32760);

    for t=1:length(t_array)
        np=t_array(t); 
        im_num=(np-1)*nslices + sl;
        
        metadata1=metadata;
        metadata1.SeriesInstanceUID = new_seuid;
        metadata1.InStackPositionNumber = sl;
        metadata1.InstanceNumber = im_num;
        metadata1.LocationsInAcquisition = nslices;
        %metadata1.LargestImagePixelValue = 6000;
        metadata1.SeriesDescription = ['RM-AW Int Cal ' metadata.SeriesDescription];
        metadata1.SeriesNumber = dataFileHeader.series.se_no*100;
        metadata1.ImagesInAcquisition = total_images;
        metadata1.ImagesInSeries = total_images;
        
        metadata1.SliceLocation = slice_demux_coord;
        metadata1.ImageLocation = slice_demux_coord;
        metadata1.ImagePositionPatient(3)= slice_demux_coord;
        
        if strcmp(dataFileHeader.image.psd_iname, 'EPI2')
            metadata1.UserData23 = num_t2;
            grad_vec = tensor_multishell(metadata.UserData24);
            if np > num_t2
                metadata1.(dicomlookup('0043','1039')) = [max_bval;8;0;0];
                diff_dir = np - num_t2;
                metadata1.UserData20=grad_vec(diff_dir,1);
                metadata1.UserData21=grad_vec(diff_dir,2);
                metadata1.UserData22=grad_vec(diff_dir,3);
            else
                metadata1.(dicomlookup('0043','1039')) = [0;8;0;0];
                metadata1.UserData20=0;
                metadata1.UserData21=0;
                metadata1.UserData22=0;
            end
        end
        
        % Setting the window and level per image
        im=img(:,:,t);
        max_im=max(im(:));
        metadata1.WindowWidth = max_im;
        metadata1.WindowCenter = max_im/2;
        metadata1.WindowValue = max_im/2;
        
        metadata1.Rows = num_rows;
        metadata1.Height =num_rows;
        metadata1.Columns=num_cols;
        metadata1.Width = num_cols;
        metadata1.PixelSpacing = pixel_spacing;
        metadata1.PixelBandwidth = pixel_bw;
        
        outfilename = sprintf('%s%05d.dcm',outname,im_num);
        %fprintf('Writing out dicom file %s', outfilename);
        dicomwrite(im, outfilename, metadata1, 'WritePrivate', true, 'CreateMode', 'copy');
        clear metadata1 outfilename;
        
        % MJM
        %movefile('/Users/matthewmiddione/Documents/MATLAB/MB_RECON/ABCD_mux_epi_recon/*.dcm','/Users/matthewmiddione/tmp/');
        
    end
end

