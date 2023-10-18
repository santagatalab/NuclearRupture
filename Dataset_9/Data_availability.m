fileinfo = dir("Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\COY_Old_LiveCell\Live_Cell_Vids_from_Alex\18_0524 Kuramochi GFP-BAF\Imsdata");
filename = {fileinfo.name};
bytelist = [fileinfo.bytes];
n = 1;
for i = 1:length(filename)
    cname = filename{i};
    if contains(cname,'.ims')
        totallist{n,1} = str2num(strtok(cname,'.'));
        totallist{n,2} = filename{i};
        totallist{n,4} = str2num(sprintf('%.0f',bytelist(i)/(1024^2)));
        n = n+1;        
    end
end

fileinfo = dir("Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\COY_Old_LiveCell\Live_Cell_Vids_from_Alex\18_0524 Kuramochi GFP-BAF\Analysis\seg");
filename = {fileinfo.name};
bytelist = [fileinfo.bytes];
n = 1;
for i = 1:length(filename)
    cname = filename{i};
    if contains(cname,'.mat')
        totallist_1{n,1} = str2num(strtok(cname,'.'));
        totallist_1{n,2} = filename{i};
        totallist_1{n,4} = str2num(sprintf('%.0f',bytelist(i)/(1024^2)));
        n = n+1;        
    end
end
