clear;
clc;

% ≠Ïload image
seq_path = 'D:\chen\(SC,DPZ,ET) & MMD\SL16-all';
dirlist = dir([seq_path '\*201*']);
output_pathname = 'D:\chen\(SC,DPZ,ET) & MMD\SL16-all';
outnum =randperm(numel(dirlist));
for directindex = 1:numel(dirlist)  % Run through every desired directory found in seq_path
    if (dirlist(directindex).isdir==1)
    input_pathname = [seq_path '\' dirlist(directindex).name];
    input_filelist = dir([input_pathname '\Chan*.tif']);
    n_frame = numel(input_filelist);
    input_filename = cell(1,n_frame);
    for inputfile_idx = 1:n_frame
        input_filename{inputfile_idx} = {input_filelist(inputfile_idx).name};
    end
    
    
    % ≠ÏSave stack
    %[token1, rem] = strtok(dirlist(directindex).name, '_'); %2015XXXX
    %[tokenmmm,rem] =strtok(rem, '_');
    %[token2, rem] = strtok(rem, '_'); %S600T6XX
    %[token3, rem] = strtok(rem, '_'); %Liao  %rem
    %[token4, rem] = strtok(rem, '_'); %SL07_ATNor01
    %stackname = [token4 tokenmmm rem '_' token2(5) token2(6) token2(7) token2(8)];
    stackname = [dirlist(directindex).name];
    for i=1:1:n_frame/2
        a_path = strcat(input_pathname, '\', input_filename{i});
        b_path = strcat(input_pathname, '\', input_filename{i+n_frame/2});
        a=imread(a_path{1});
        b=imread(b_path{1});
        thg = cat(3,b,zeros(size(b)),b);
        shg = cat(3,zeros(size(a)),a,zeros(size(a)));
        img = thg + shg;
        imwrite(img,[output_pathname,'\',num2str(outnum(directindex)),'.tif'],'WriteMode','append','Compression','none');
    end
   
    end
end
