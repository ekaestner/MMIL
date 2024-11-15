
ttt = niftiread('/home/darienzo/Documents/MATLAB/spm12/toolbox/cat12/templates_MNI152NLin2009cAsym/aal3.nii');
slc_130 = rot90(squeeze(ttt(:,130,:)));
figure(); imagesc(slc_130)
tabulate(slc_130(:))

tt3 = interp3(ttt,'method','nearest')