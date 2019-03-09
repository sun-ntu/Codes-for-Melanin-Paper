

%% ----- do statistics -----

%% Load Image
% n_TPF_file   %read the name of TPF file
% n_fTHG_file  %read the name of fTHG file
% n_bTHG_file  %read the name of bTHG file
% n_BF_file    %read the name of BF file
%
% img_TPF : read file as img fmt (512x512)
% img_TPF : read file as img fmt (512x512)
% img_fTHG : read file as img fmt (512x512)
% img_bTHG : read file as img fmt (512x512)
% img_BF : read file as img fmt (512x512)


clc
close all
clear all


%Input number of file to select data

no_HG = input('Input the number of HG images:  ', 's');
no_BF = input('Input the number of BF images:  ', 's');

n_TPF_file0 = 'Image'; 
n_TPF_file1 = no_HG; 
n_TPF_file2 = '_ch01'; 
n_TPF_file3 = '.tif';
n_TPF_file = [n_TPF_file0 n_TPF_file1 n_TPF_file2 n_TPF_file3];


n_fTHG_file0 = 'Image'; 
n_fTHG_file1 = no_HG; 
n_fTHG_file2 = '_ch00'; 
n_fTHG_file3 = '.tif';
n_fTHG_file = [n_fTHG_file0 n_fTHG_file1 n_fTHG_file2 n_fTHG_file3];


n_bTHG_file0 = 'Image'; 
n_bTHG_file1 = no_HG; 
n_bTHG_file2 = '_ch02'; 
n_bTHG_file3 = '.tif';
n_bTHG_file = [n_bTHG_file0 n_bTHG_file1 n_bTHG_file2 n_bTHG_file3];


n_BF_file0 = 'Image'; 
n_BF_file1 = no_BF; 
n_BF_file2 = '_ch00'; 
n_BF_file3 = '.tif';
n_BF_file = [n_BF_file0 n_BF_file1 n_BF_file2 n_BF_file3];



% Use mouse to select file 

% n_TPF_file  = uigetfile('*.tif'); %read TPF image
% n_fTHG_file = uigetfile('*.tif'); %read fTHG image
% n_bTHG_file = uigetfile('*.tif'); %read bTHG image
% n_BF_file   = uigetfile('*.tif'); %read BF image

display(n_TPF_file);
display(n_fTHG_file);
display(n_bTHG_file);
display(n_BF_file);

img_TPF = imread(n_TPF_file) ;       % 載入一張螢光影像 (TPF)
img_fTHG = imread(n_fTHG_file) ;     % 載入一張正向三倍頻影像 (fTHG)
img_bTHG = imread(n_bTHG_file) ;     % 載入一張背向三倍頻影像 (bTHG)
img_BF = imread(n_BF_file) ;         % 載入一張BF影像

%% Set 'fn_Data' :  name of file 
% 'fn_Data' : record the name of data
%fn = [fn_Data  fn_Time fn_Content]


no_LengthOfName=8; %ex. 'Image251_ch01.tif' => 'Image251'

fn_Data = n_fTHG_file(1:no_LengthOfName);

%% Set 'fn_Time' : Time of File name(1)
% ex. fn_Time = _20130723_1636

%edit: 20130808  filename: 'fn_Data' _'fn_Time' _'fn_Content' ...
fn_Time = clock; 
fn_Time1= '_';

if fn_Time(2) < 10 
fn_Time2 = '0';
else
fn_Time2 = '';      
end

if fn_Time(3) < 10 
fn_Time3 = '0';
else
fn_Time3 = '';      
end


if fn_Time(4) < 10 
fn_Time4 = '0';
else
fn_Time4 = '';      
end

if fn_Time(5) < 10 
fn_Time5 = '0';
else
fn_Time5 = '';      
end

fn_Time = [fn_Time1 num2str(fn_Time (1,1)) fn_Time2 num2str(fn_Time (1,2)) fn_Time3 ...
           num2str(fn_Time (1,3)) fn_Time1 fn_Time4 num2str(fn_Time (1,4)) fn_Time5 ...
           num2str(fn_Time (1,5))];
       

%% Select ROI
% roi_region: Roi region
% Mask_Show_img: selected region


fTHG_Show_img = zeros(512,512,3);
fTHG_Show_img = uint16(fTHG_Show_img);
fTHG_Show_img(:,:,1) = img_fTHG*16;
fTHG_Show_img(:,:,2) = 0;
fTHG_Show_img(:,:,3) = img_fTHG*16;


bTHG_Show_img = zeros(512,512,3);
bTHG_Show_img = uint16(bTHG_Show_img);
bTHG_Show_img(:,:,1) = img_bTHG*16;
bTHG_Show_img(:,:,2) = 0;
bTHG_Show_img(:,:,3) = img_bTHG*16;



TPF_Show_img = zeros(512,512,3);
TPF_Show_img = uint16(TPF_Show_img);
TPF_Show_img(:,:,1) = img_TPF*16;
TPF_Show_img(:,:,2) = img_TPF*16;
TPF_Show_img(:,:,3) = 0;




img_TPF = imread(n_TPF_file) ;       % 載入一張螢光影像 (TPF)


roi_region = logical(zeros(512,512));

BF_Show_img = zeros(512,512,3);
BF_Show_img = uint16(BF_Show_img);
BF_Show_img(:,:,1) = img_BF*16;
BF_Show_img(:,:,2) = img_BF*16;
BF_Show_img(:,:,3) = img_BF*16;

Mask_Show_img = zeros(512,512,3);
Mask_Show_img = uint16(Mask_Show_img);



% imshow(3*fTHG_Show_img);
Flag = 'A';
%Flag_show = input('What kind of image? (B/T)(B:bf+fTHG, T:fTHG):  ', 's');

while Flag == 'A'|'D'
    
    Flag_show = input('What kind of image? (B/T/return:quit)(B:bf+fTHG, T:fTHG):  ', 's');
    
    switch Flag
        
        case 'A'
            
            switch Flag_show
                case 'B'
                    
                    imshow(BF_Show_img + Mask_Show_img);
                    
                    roi_region = roi_region|roipoly;
                    
                    Mask_Show_img(:,:,1) = roi_region*2.^15;
                    imshow(BF_Show_img + Mask_Show_img);
                case 'T'
                    
                    imshow(3*fTHG_Show_img + Mask_Show_img);
                    
                    roi_region = roi_region|roipoly;
                    
                    Mask_Show_img(:,:,1) = roi_region*2.^15;
                    imshow(3*fTHG_Show_img + Mask_Show_img);
            end
            
            
        case 'D'
            
            switch Flag_show
                case 'B'
                    
                    imshow(BF_Show_img + Mask_Show_img);
                    
                    roi_region = roi_region&(~(roi_region & roipoly));
                    
                    Mask_Show_img(:,:,1) = roi_region*2.^15;
                    imshow(BF_Show_img + Mask_Show_img);
                    
                case 'T'
                    
                    imshow(3*fTHG_Show_img + Mask_Show_img);
                    
                    
                    roi_region = roi_region&(~(roi_region & roipoly));
                    
                    Mask_Show_img(:,:,1) = roi_region*2.^15;
                    imshow(3*fTHG_Show_img + Mask_Show_img);
            end
            
            
    end
    
    
    
    
    Flag = input('Do you want to continue? ([A]:add / [D]:delete / Return:quit) ', 's');
    
end


%% Import Data
%Read data, Set value to
% fTHG_TPF_dot
% bTHG_TPF_dot
% fTHG_THG_dot
% bTHG_THG_dot
% Sat_fTHG_TPF_dot
% Sat_fTHG_THG_dot
% Sat_bTHG_TPF_dot
% Sat_bTHG_THG_dot


%1) ROI points
% fTHG
% TPF:  fTHG_TPF_dot
% THG:  fTHG_THG_dot
% bTHG
% TPF:  bTHG_TPF_dot
% THG:  bTHG_THG_dot

%2)saturation point
% fTHG
% TPF:  Sat_fTHG_TPF_dot
% THG:  Sat_fTHG_THG_dot
% bTHG
% TPF:  Sat_bTHG_TPF_dot
% THG:  Sat_bTHG_THG_dot

roi_pixel_count = 0; % count number of ROI pixel


selected_roi_region_idx = find ( roi_region == 1 );
no_roi_region = size ( selected_roi_region_idx );

fTHG_TPF_dot = zeros(no_roi_region);
fTHG_THG_dot = zeros(no_roi_region);

bTHG_TPF_dot = zeros(no_roi_region);
bTHG_THG_dot = zeros(no_roi_region);

Sat_fTHG_count=0;
Sat_bTHG_count=0;

% Calculate the number of saturation fTHG, bTHG pixels

for i=1:512
    for j=1:512
        if roi_region(i,j) ~= 0
            
            if img_fTHG(i,j) == 4095
                Sat_fTHG_count = Sat_fTHG_count+1;
            end
            
            if img_bTHG(i,j) == 4095
                Sat_bTHG_count = Sat_bTHG_count+1;
            end
            
        end
    end
end


display(Sat_fTHG_count);
display(Sat_bTHG_count);

Sat_fTHG_TPF_dot = zeros (Sat_fTHG_count);
Sat_fTHG_THG_dot = zeros (Sat_fTHG_count);
Sat_bTHG_TPF_dot = zeros (Sat_bTHG_count);
Sat_bTHG_THG_dot = zeros (Sat_bTHG_count);


% set value for
% fTHG_TPF_dot
% bTHG_TPF_dot
% fTHG_THG_dot
% bTHG_THG_dot
% Sat_fTHG_TPF_dot
% Sat_fTHG_THG_dot
% Sat_bTHG_TPF_dot
% Sat_bTHG_THG_dot



for i=1:512
    for j=1:512
        if roi_region(i,j) ~= 0
            if img_fTHG(i,j)==4095
                Sat_fTHG_TPF_dot(Sat_fTHG_count+1) = img_TPF(i,j); %sPIC 存所有擁有飽和THG的TPF的值
                Sat_fTHG_THG_dot(Sat_fTHG_count+1) = img_fTHG(i,j);
            end
            
            if img_bTHG(i,j)==4095
                Sat_bTHG_TPF_dot(Sat_bTHG_count+1) = img_TPF(i,j); %sPIC 存所有擁有飽和THG的TPF的值
                Sat_bTHG_THG_dot(Sat_bTHG_count+1) = img_fTHG(i,j);
            end
            
            
            fTHG_TPF_dot(roi_pixel_count+1) = img_TPF(i,j);%TPF
            bTHG_TPF_dot(roi_pixel_count+1) = img_TPF(i,j);%TPF
            fTHG_THG_dot(roi_pixel_count+1) = img_fTHG(i,j);%fTHG
            bTHG_THG_dot(roi_pixel_count+1) = img_bTHG(i,j);%bTHG
            
            roi_pixel_count = roi_pixel_count+1;
            
        end
    end
end


%% Find background THG
% %將所有TPF為0時的THG平均起來
% %bg_fTHG : mean of forward  THG pixel value of data with their TPF = 0
% %bg_bTHG : mean of backward THG pixel value of data with their TPF = 0
%
% no_zeroTPF = find(fTHG_TPF_dot(:) == 0);
%
% bg_fTHG = mean(fTHG_THG_dot(no_zeroTPF));
% bg_bTHG = mean(bTHG_THG_dot(no_zeroTPF));
%
%

%% Disperse into partitions
%TPF>=1 開始
%Psep: separation distance of TPF pixel value
%fTHG_THG_mean: averaged forward  THG value in each partition
%bTHG_THG_dotY: averaged backward THG value in each partition

%Sepa_TPF_dot : structure saving indexs of TPF data in each partition


%fTHG_ST : standard deviation of fTHG value in each partition
%bTHG_ST : standard deviation of bTHG value in each partition

% Psep = input('Enter a Number for number of pixel in a seperation : ' ); % 輸入一數值將螢光強度分群，其輸入數字為二的整數倍。 Ex:2、4、6....

Psep = 50;


k = 1 ;

for I = 1:Psep:4095
    
    [index_TPF,y] = find((fTHG_TPF_dot(:) >= I)&(fTHG_TPF_dot(:) <I+Psep)) ; % [x,y]為那些值的[row,column]的值,y已經被固定了，故不care。
    Sepa_TPF_dot(k) = struct('X_Axis' , index_TPF) ;
    k = k+1 ;
    
end


if mod(4094,Psep)==0
    Nsep= (4094-mod(4094,Psep))/Psep;
else
    Nsep=(4094-mod(4094,Psep))/Psep+1;
end



for r = 1:Nsep
    
    TPF_mean(r) = mean(fTHG_TPF_dot(getfield(Sepa_TPF_dot(1,r) , 'X_Axis'))) ; % 求各分群中其螢光強度的平均值
    
end

for t = 1:Nsep
    
    fTHG_THG_mean(t) = mean (fTHG_THG_dot(getfield(Sepa_TPF_dot(1,t) , 'X_Axis'))) ; % 求相對各分群中螢光強度的有效三倍頻強度(平均值)
    bTHG_THG_mean(t) = mean (bTHG_THG_dot(getfield(Sepa_TPF_dot(1,t) , 'X_Axis'))) ; % 求相對各分群中螢光強度的有效三倍頻強度(平均值)
    fTHG_ST(t) = std (fTHG_THG_dot(getfield(Sepa_TPF_dot(1,t) , 'X_Axis')));
    bTHG_ST(t) = std (bTHG_THG_dot(getfield(Sepa_TPF_dot(1,t) , 'X_Axis')));
    
end


%% Calculate contribution of saturation pixels (fTHG)
%Sepa_Sat_fTHG_TPF_dot : structure saving indexs of Sat data in each TPF partition
%SumSat_fTHG: 每個分隔飽和值pixel值加總
%SumAll_fTHG: 每個分隔內全部pixel值加總
%SumR_fTHG  : 每個分隔內飽和值所貢獻pixel值比例(in%)



k=1 ;

for I = 1:Psep:4095
    
    [index_Sat_fTHG_TPF_dot,y] = find((Sat_fTHG_TPF_dot(:) >= I)&(Sat_fTHG_TPF_dot(:) <I+Psep)) ; % [x,y]為那些值的[row,column]的值,y已經被固定了，故不care。
    Sepa_Sat_fTHG_TPF_dot(k) = struct('X_Axis' , index_Sat_fTHG_TPF_dot) ;
    k = k+1 ;
end


for r = 1:Nsep
    
    Sum_Sat_fTHG_dot(r) = sum(Sat_fTHG_THG_dot(getfield(Sepa_Sat_fTHG_TPF_dot(1,r) , 'X_Axis')));
    Sum_All_fTHG_dot(r) = sum(fTHG_THG_dot(getfield(Sepa_TPF_dot(1,r) , 'X_Axis')));
   
    
end

 Sat_Rate_fTHG = Sum_Sat_fTHG_dot./Sum_All_fTHG_dot*100;


%依據飽和比例分群，某比例以上呈現紅色，以下呈現綠色

no_red_SumR_fTHG = 0;
no_green_SumR_fTHG=0;



for r = 1:Nsep
    
    if Sat_Rate_fTHG(r)>20
        
        red_Sat_Rate_fTHG(no_red_SumR_fTHG+1) = Sat_Rate_fTHG(r);
        red_Sat_Rate_fTHG(no_red_SumR_fTHG+1) = TPF_mean(r);
        no_red_SumR_fTHG = no_red_SumR_fTHG+1;
        
    else
        
        green_Sat_Rate_fTHG(no_green_SumR_fTHG+1) = Sat_Rate_fTHG(r);
        green_Sat_Rate_fTHG(no_green_SumR_fTHG+1) = TPF_mean(r);
        no_green_SumR_fTHG = no_green_SumR_fTHG+1;
        
    end
    
end

%% Calculate contribution of saturation pixels (bTHG)
%Sepa_Sat_bTHG_TPF_dot : structure saving indexs of Sat data in each TPF partition
%SumSat_bTHG: 每個分隔飽和值pixel值加總
%SumAll_bTHG: 每個分隔內全部pixel值加總
%SumR_bTHG  : 每個分隔內飽和值所貢獻pixel值比例(in%)



k=1 ;

for I = 1:Psep:4095
    
    [index_Sat_bTHG_TPF_dot,y] = find((Sat_bTHG_TPF_dot(:) >= I)&(Sat_bTHG_TPF_dot(:) <I+Psep)) ; % [x,y]為那些值的[row,column]的值,y已經被固定了，故不care。
    Sepa_Sat_bTHG_TPF_dot(k) = struct('X_Axis' , index_Sat_bTHG_TPF_dot) ;
    k = k+1 ;
    
end

for r = 1:Nsep
    
    Sum_Sat_bTHG_dot(r) = sum(Sat_bTHG_THG_dot(getfield(Sepa_Sat_bTHG_TPF_dot(1,r) , 'X_Axis')));
    Sum_All_bTHG_dot(r) = sum(fTHG_THG_dot(getfield(Sepa_TPF_dot(1,r) , 'X_Axis')));
    
end

Sat_Rate_bTHG = Sum_Sat_bTHG_dot./Sum_All_bTHG_dot*100;

%依據飽和比例分群，某比例以上呈現紅色，以下呈現綠色


no_red_SumR_bTHG = 0;
no_green_SumR_bTHG = 0;



for r = 1:Nsep
    if Sat_Rate_bTHG(r)>20
        red_Sat_Rate_bTHG(no_red_SumR_bTHG+1) = Sat_Rate_bTHG(r);
        red_Sat_Rate_bTHG(no_red_SumR_bTHG+1) = TPF_mean(r);
        no_red_SumR_bTHG = no_red_SumR_bTHG+1;
    else
        green_Sat_Rate_bTHG(no_green_SumR_bTHG+1) = Sat_Rate_bTHG(r);
        green_Sat_Rate_bTHG(no_green_SumR_bTHG+1) = TPF_mean(r);
        no_green_SumR_bTHG = no_green_SumR_bTHG+1;
    end
    
end



%% ------ present -------

%% Plot dot figure (fTHG)
%plot "fTHG pixel vs TPF pixel dot" figure

%Set XY parameters

%Hold 

%Plot
figure,
plot(fTHG_TPF_dot(:),fTHG_THG_dot(:),'.b');


%Set title 
% Title = 'Empirical Equation';
% size_Title = 15 %set
% T = title(Title,'Fontsize',size_Title);
% 
% MA    = 'Filename: ';
% RW    = 'Raw Data';
% N_fTHG= 'fTHG';
% Title = [MA fn_Data 10 RW 10 N_fTHG];
% title(Title)

%Set legend 
size_legend = 19; %set


legend_Data = 'Data Point';
L = legend(legend_Data);
set(L,'FontSize',size_legend);%set


%Set ticks
size_ticks = 20; 
set(gca,'FontSize',size_ticks) %Set font size to tick lables, the same as legend 

%Set label
size_label = 24; 
xlabel('TPF (pixel value) ','Fontsize',size_label)
ylabel('THG (pixel value)','Fontsize',size_label)

%box
box on

%grid

%Set axis
xlim([0 4200])
ylim([0 4200])

%Save 







%% Plot dot figure (bTHG)
%plot "bTHG pixel value vs TPF pixel value dot" figure

%Set XY parameters

%Hold 

%Plot
figure,
plot(bTHG_TPF_dot(:),bTHG_THG_dot(:),'.m');


%Set title 
% title('Melanin Analysis (mean value)')
% MA    = 'Filename: ';
% RW    = 'Raw Data';
% N_bTHG= 'bTHG';
% Title = [MA fn_Data 10 RW 10 N_bTHG];
% title(Title)


%Set legend 
size_legend = 19; %set


legend_Data = 'Data Point';
L = legend(legend_Data);
set(L,'FontSize',size_legend);%set


%Set ticks
size_ticks = 20; 
set(gca,'FontSize',size_ticks) %Set font size to tick lables, the same as legend 

%Set label
size_label = 24; 
xlabel('TPF (pixel value)','Fontsize',size_label)
ylabel('THG (pixel value)','Fontsize',size_label)

%box
box on

%grid

%Set axis
xlim([0 4200])
ylim([0 4200])


%Save 


%%  Plot error+mean figure (fTHG)
%fTHG
%plot data point with error bar

%Set XY parameters

N_TPF_mean=size(TPF_mean);
N_TPF_mean=N_TPF_mean(2);

%Set data range for fitting

%Plot
figure,
errorbar(TPF_mean,fTHG_THG_mean,fTHG_ST,'bo')  % 畫出標準差


% %Set title 
% MA    = 'Filename: ';
% AvgErr='Average with error';
% SP0   = 'Number of TPF Pixels per seperation = ';
% SP1   = [SP0 num2str(Psep) ];
% 
% Title = [MA fn_Data 10 AvgErr 10 SP1 10 N_fTHG];
% title(Title)

%Set legend 
size_legend = 19; %set

legend_Data = 'Mean Value of Data in Each Partition';
L = legend(legend_Data);
set(L,'FontSize',size_legend);%set


%Set ticks
size_ticks = 20; 
set(gca,'FontSize',size_ticks) %Set font size to tick lables, the same as legend 

%Set label
size_label = 24; 
xlabel('TPF (pixel value)','Fontsize',size_label)
ylabel('THG (pixel value)','Fontsize',size_label)

%box
box on

%grid

%Set axis
xlim([0 4200])
ylim([0 4200])


%%  Plot error+mean figure (bTHG)
%bTHG
%plot data point with error bar

%Set XY parameters

N_TPF_mean=size(TPF_mean);
N_TPF_mean=N_TPF_mean(2);

%Set data range for fitting

%Plot
figure,
errorbar(TPF_mean,bTHG_THG_mean,bTHG_ST,'mo')   % 畫出標準差


%Set title 
% MA    = 'Filename: ';
% AvgErr='Average with error';
% SP0   = 'Number of TPF Pixels per seperation = ';
% SP1   = [SP0 num2str(Psep) ];
% Title = [MA fn_Data 10 AvgErr 10 SP1 10 N_bTHG];
% title(Title)

%Set legend 
size_legend = 19; %set

legend_Data = 'Mean Value of Data in Each Partition';
L = legend(legend_Data);
set(L,'FontSize',size_legend);%set


%Set ticks
size_ticks = 20; 
set(gca,'FontSize',size_ticks) %Set font size to tick lables, the same as legend 

%Set label
size_label = 24; 
xlabel('TPF (pixel value)','Fontsize',size_label)
ylabel('THG (pixel value)','Fontsize',size_label)

%box
box on

%grid

%Set axis
xlim([0 4200])
ylim([0 4200])





%% Plot Saturation Contribution  (fTHG)



if no_green_SumR_fTHG ~= 0
    
    N_green_Sat_Rate_fTHG=size(green_Sat_Rate_fTHG);
    N_green_Sat_Rate_fTHG=N_green_Sat_Rate_fTHG(2);
    
    
end


if no_red_SumR_fTHG ~= 0
    
    %Pixel2Concentration
    N_R_SatX=size(red_Sat_Rate_fTHG);
    N_R_SatX=N_R_SatX(2);
    
    
end




%     Show above, below 20% figure,

%     if no_green_SumR_fTHG ~=0
%         plot(green_Sat_Rate_fTHG,green_Sat_Rate_fTHG,'og')
%     end
%     hold on
%
%     if no_red_SumR_fTHG ~=0
%         plot(red_Sat_Rate_fTHG,red_Sat_Rate_fTHG,'or')
%     end
%     hold off
%



%Plot
figure,
plot(TPF_mean,Sat_Rate_fTHG,'b*');


%Set title 
% MA    = 'Filename: ';
% SP0   = 'Number of TPF Pixels per seperation = ';
% SP1   = [SP0 num2str(Psep) ];
% Title0=' Saturation data pixel%  ';
% Title = [MA fn_Data 10 Title0 10 SP1 10];
% title(Title)

%Set legend 
size_legend = 19; %set
legend_Data = 'Mean Saturation Rate';
L = legend(legend_Data);
set(L,'FontSize',size_legend);%set


%Set ticks
size_ticks = 20; 
set(gca,'FontSize',size_ticks) %Set font size to tick lables, the same as legend 

%Set label
size_label = 24; 
xlabel('TPF (pixel value)','Fontsize',size_label)
ylabel('Saturation P%','Fontsize',size_label)

%box
box on

%grid

%Set axis
xlim([0 4200])





%% Plot Saturation Contribution  (bTHG)



    if no_green_SumR_bTHG ~= 0

        N_green_Sat_Rate_bTHG=size(green_Sat_Rate_bTHG);
        N_green_Sat_Rate_bTHG=N_green_Sat_Rate_bTHG(2);


    end


    if no_red_SumR_bTHG ~= 0

        %Pixel2Concentration
        N_R_SatX=size(red_Sat_Rate_bTHG);
        N_R_SatX=N_R_SatX(2);

    end



%     Show above, below 20% figure,

%     figure,
%     if N_green_Sat_Rate_bTHG ~=0
%         plot(green_Sat_Rate_bTHG,green_Sat_Rate_bTHG,'og')
%     end
%     hold on
%
%     if N_red_Sat_Rate_bTHG ~=0
%         plot(red_Sat_Rate_bTHG,red_Sat_Rate_bTHG,'or')
%     end
%     hold off



%Plot
figure,
plot(TPF_mean,Sat_Rate_bTHG,'m*');


%Set title 
% MA    = 'Filename: ';
% SP0   = 'Number of TPF Pixels per seperation = ';
% SP1   = [SP0 num2str(Psep) ];
% Title0=' Saturation data pixel%  ';
% Title = [MA fn_Data 10 Title0 10 SP1 10];
% title(Title)

%Set legend 
size_legend = 19; %set
legend_Data = 'Mean Saturation Rate';
L = legend(legend_Data);
set(L,'FontSize',size_legend);%set


%Set ticks
size_ticks = 20; 
set(gca,'FontSize',size_ticks) %Set font size to tick lables, the same as legend 

%Set label
size_label = 24; 
xlabel('TPF (pixel value)','Fontsize',size_label)
ylabel('Saturation P%','Fontsize',size_label)

%box
box on

%grid

%Set axis
xlim([0 4200])
