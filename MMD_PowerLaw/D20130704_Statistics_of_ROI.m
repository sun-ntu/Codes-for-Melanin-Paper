
%% 
% D20130704_Statistics_of_ROI_20190308.m �����������F�����n����X �Τ����n������ (for submitting) ��X�G���4
% D20130704_Statistics_of_ROI_20190118.m �������ʨ}�ץ��F��4���B��                               ��X�G���123456
% D20130704_Statistics_of_ROI.m          ���������¥����Ѫ���                                    ��X�G���123456

%% Load Data

clc
clear all
close all

N_file='D20130704_m20130807ROI_20130808_1154_SumROI_Selected_Variable';
load D20130704_m20130807ROI_20130808_1154_SumROI_Selected_Variable;

%% Set 'fn_Data' :  name of file 
% 'fn_Data' : record the name of data
%fn = [fn_Data  fn_Time fn_Content]

fn_Data = N_file;


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
       


%% Disperse into partitions
%TPF>=1 �}�l
%Psep: separation distance of TPF pixel value
%fTHG_THG_mean: averaged forward  THG value in each partition
%bTHG_THG_dotY: averaged backward THG value in each partition

%Sepa_TPF_dot : structure saving indexs of TPF data in each partition


%fTHG_ST : standard deviation of fTHG value in each partition
%bTHG_ST : standard deviation of bTHG value in each partition

% Psep = input('Enter a Number for number of pixel in a seperation : ' ); % ��J�@�ƭȱN�å��j�פ��s�A���J�Ʀr���G����ƭ��C Ex:2�B4�B6....

Psep = 50;


k = 1 ;

for I = 1:Psep:4095
    
    [index_TPF,y] = find((fTHG_TPF_dot(:) >= I)&(fTHG_TPF_dot(:) <I+Psep)) ; % [x,y]�����ǭȪ�[row,column]����,y�w�g�Q�T�w�F�A�G��care�C
    Sepa_TPF_dot(k) = struct('X_Axis' , index_TPF) ;
    k = k+1 ;
    
end


if mod(4094,Psep)==0
    Nsep= (4094-mod(4094,Psep))/Psep;
else
    Nsep=(4094-mod(4094,Psep))/Psep+1;
end



for r = 1:Nsep
    
    TPF_mean(r) = mean(fTHG_TPF_dot(getfield(Sepa_TPF_dot(1,r) , 'X_Axis'))) ; % �D�U���s����å��j�ת�������
    
end

for t = 1:Nsep
    
    fTHG_THG_mean(t) = mean (fTHG_THG_dot(getfield(Sepa_TPF_dot(1,t) , 'X_Axis'))) ; % �D�۹�U���s���å��j�ת����ĤT���W�j��(������)
    bTHG_THG_mean(t) = mean (bTHG_THG_dot(getfield(Sepa_TPF_dot(1,t) , 'X_Axis'))) ; % �D�۹�U���s���å��j�ת����ĤT���W�j��(������)
    fTHG_ST(t) = std (fTHG_THG_dot(getfield(Sepa_TPF_dot(1,t) , 'X_Axis')));
    bTHG_ST(t) = std (bTHG_THG_dot(getfield(Sepa_TPF_dot(1,t) , 'X_Axis')));
    bTHG_SE(t) = bTHG_ST(t)/sqrt(length(Sepa_TPF_dot(1,t).X_Axis));                  % modified by Ming-Liang 20190111 from Standard Deviation to Standard Error
end


%% Set values to Saturation data

Sat_fTHG_count = 0;
Sat_bTHG_count = 0;


for i= 1: Tot_no_Pixel
    
    if fTHG_THG_dot(i) == 4095
        
        Sat_fTHG_TPF_dot(Sat_fTHG_count+1) = fTHG_TPF_dot(i); %sPIC �s�Ҧ��֦����MTHG��TPF����
        Sat_fTHG_THG_dot(Sat_fTHG_count+1) = fTHG_THG_dot(i);
        Sat_fTHG_count = Sat_fTHG_count + 1;
        
    end
    
    if bTHG_THG_dot(i) == 4095
        
        Sat_bTHG_TPF_dot(Sat_bTHG_count+1) = bTHG_TPF_dot(i); %sPIC �s�Ҧ��֦����MTHG��TPF����
        Sat_bTHG_THG_dot(Sat_bTHG_count+1) = bTHG_THG_dot(i);
        Sat_bTHG_count = Sat_bTHG_count + 1;
   
    end
    
end


display (Sat_fTHG_count);
display (Sat_bTHG_count);


%% Calculate contribution of saturation pixels (fTHG)
%Sepa_Sat_fTHG_TPF_dot : structure saving indexs of Sat data in each TPF partition
%SumSat_fTHG: �C�Ӥ��j���M��pixel�ȥ[�`
%SumAll_fTHG: �C�Ӥ��j������pixel�ȥ[�`
%SumR_fTHG  : �C�Ӥ��j�����M�ȩҰ^�mpixel�Ȥ��(in%)



k=1 ;

for I = 1:Psep:4095
    
    [index_Sat_fTHG_TPF_dot,y] = find((Sat_fTHG_TPF_dot(:) >= I)&(Sat_fTHG_TPF_dot(:) <I+Psep)) ; % [x,y]�����ǭȪ�[row,column]����,y�w�g�Q�T�w�F�A�G��care�C
    Sepa_Sat_fTHG_TPF_dot(k) = struct('X_Axis' , index_Sat_fTHG_TPF_dot) ;
    k = k+1 ;
end


for r = 1:Nsep
    
    Sum_Sat_fTHG_dot(r) = sum(Sat_fTHG_THG_dot(getfield(Sepa_Sat_fTHG_TPF_dot(1,r) , 'X_Axis')));
    Sum_All_fTHG_dot(r) = sum(fTHG_THG_dot(getfield(Sepa_TPF_dot(1,r) , 'X_Axis')));
   
    
end

 Sat_Rate_fTHG = Sum_Sat_fTHG_dot./Sum_All_fTHG_dot*100;


%�̾ڹ��M��Ҥ��s�A�Y��ҥH�W�e�{����A�H�U�e�{���

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
%SumSat_bTHG: �C�Ӥ��j���M��pixel�ȥ[�`
%SumAll_bTHG: �C�Ӥ��j������pixel�ȥ[�`
%SumR_bTHG  : �C�Ӥ��j�����M�ȩҰ^�mpixel�Ȥ��(in%)



k=1 ;

for I = 1:Psep:4095
    
    [index_Sat_bTHG_TPF_dot,y] = find((Sat_bTHG_TPF_dot(:) >= I)&(Sat_bTHG_TPF_dot(:) <I+Psep)) ; % [x,y]�����ǭȪ�[row,column]����,y�w�g�Q�T�w�F�A�G��care�C
    Sepa_Sat_bTHG_TPF_dot(k) = struct('X_Axis' , index_Sat_bTHG_TPF_dot) ;
    k = k+1 ;
    
end

for r = 1:Nsep
    
    Sum_Sat_bTHG_dot(r) = sum(Sat_bTHG_THG_dot(getfield(Sepa_Sat_bTHG_TPF_dot(1,r) , 'X_Axis')));
    Sum_All_bTHG_dot(r) = sum(fTHG_THG_dot(getfield(Sepa_TPF_dot(1,r) , 'X_Axis')));
    
end

Sat_Rate_bTHG = Sum_Sat_bTHG_dot./Sum_All_bTHG_dot*100;

%�̾ڹ��M��Ҥ��s�A�Y��ҥH�W�e�{����A�H�U�e�{���


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

%%  Plot error+mean figure (bTHG)
%bTHG
%plot data point with error bar

%Set XY parameters

N_TPF_mean=size(TPF_mean);
N_TPF_mean=N_TPF_mean(2);

%Set data range for fitting

%Plot
figure,
% errorbar(TPF_mean,bTHG_THG_mean,bTHG_ST,'mo')   % �e�X�зǮt
errorbar(TPF_mean,bTHG_THG_mean,bTHG_SE,'mo')   % modified by Ming-Liang 20190111 from Standard Deviation to Standard Error
title('THG-TPF')

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
