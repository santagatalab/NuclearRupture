clear
basedir = 'Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\COY_Old_LiveCell\Live_Cell_Vids_from_Alex\18_0524 Kuramochi GFP-BAF';
savedir = [basedir '\Analysis\gfp'];
mkdir(savedir)
%Imaris Reading main
for index = 4
    fileObj = ImarisReader([basedir '\Imsdata\18_054 Kuramochi GFP-BAF_' num2str(index) '_18_054 Kuramochi GFP-BAF.ims']);
    
    %gfp focus
    cIdx = 1;
    tIdx = 1;
    image_gfplist = fileObj.DataSet.GetDataVolume(cIdx, tIdx);
    [mx,my,mz] = size(image_gfplist);
    focuslist = [];
    for z = 1:mz
        fm = fmeasure(image_gfplist(:,:,z),'SFIL');
        focuslist = [focuslist fm];
    end
    [~,focusz] = max(focuslist);
    for tIdx = 1:fileObj.DataSet.SizeT-1
        image_gfpc = fileObj.DataSet.GetDataVolume(cIdx, tIdx);
        image_gfp_allt(:,:,tIdx) = image_gfpc(:,:,focusz).';
    end
    options.overwrite = true;
    saveastiff(image_gfp_allt,[savedir '\' num2str(index) '_gfp.tiff'],options)
end