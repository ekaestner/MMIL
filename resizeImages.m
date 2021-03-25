height=400;
width=400;

mainDirectory = 'C:\\Users\\iasp\\Downloads\\OneDrive_2020-10-23\\all stimuli\\Set_%s\\';
folder = ['C' 'D' 'E' 'F'];

for i=1:length(folder)
    stimuliFolder = sprintf(mainDirectory,folder(i));
    filename = fullfile(stimuliFolder, '*.jpg');
    allFiles = dir(filename);
    
    for j=1:length(allFiles)
        imageFilename = fullfile(allFiles(j).folder, allFiles(j).name);
        [filepath,name,ext] = fileparts(imageFilename);
        image = imread(imageFilename);
        
        if isequal(size(image), [400 400 3])
            movefile(imageFilename, strcat(stimuliFolder, name, '_400x400', ext));
        else
           imageResized = imresize(image, [height, width]);
           imwrite(imageResized, strcat(stimuliFolder, name, '_400x400', ext));
        end
    end
end