%% 
% D20130704_PolyRegression_20131028_20190308.m 此版本移除了不必要的輸出 及不必要的註解 (for submitting) 輸出：原圖37
% D20130704_PolyRegression_20131028.m          此版本為威民提供版本                                    輸出：原圖123456789

%% ==== THG ratio (y) vs. MMD (x) ====
   

%% ----- do statistics -----

clc
clear all
close all

%% Set parameters 

%
size_FitOne = 16;
%
bg_THGratio = 1.06;
%
zero_MMD = '20131028';
zero_MMD_TPF = 390.757; %20131028
MMD_Calibration = '20131028';


%% Load data

%data source: m_20137191925_Statistics_D20130704 

%bTHG
%plot data point with error bar


load D20130704_m20130807ROI_20130808_1154_SumROI_Selected_Variable_20130808_1611_statistics;

%% Set 'fn_Data' :  name of file 
% 'fn_Data' : record the name of data
%fn = [fn_Data  fn_Time fn_Content]



fn_Data = 'D20130704_m20130807ROI_20130808_1154_SumROI_Selected_Variable_20130808_1611_statistics';


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

       

       

%% Disperse into partitions (THG vs. MMD)
%20130805
%TPF>=1 開始
%Psep: separation distance of TPF pixel value
%fTHG_THG_mean: averaged forward  THG value in each partition
%bTHG_THG_dotY: averaged backward THG value in each partition

%Sepa_TPF_dot : structure saving indexs of TPF data in each partition


%fTHG_ST : standard deviation of fTHG value in each partition
%bTHG_ST : standard deviation of bTHG value in each partition

% Psep = input('Enter a Number for number of pixel in a seperation : ' ); % 輸入一數值將螢光強度分群，其輸入數字為二的整數倍。 Ex:2、4、6....




[index_TPF,y] = find((fTHG_TPF_dot(:) < zero_MMD_TPF)) ; 

Sepa_TPF_dot(1) = struct('X_Axis' , index_TPF) ;
 
Psep = 50;
k = 2 ;




if mod(4095-zero_MMD_TPF,Psep)==0
    Nsep = (4095-zero_MMD_TPF)/Psep + 1; % + 1 for background 
else
    Nsep=(4095-zero_MMD_TPF -mod(4095-zero_MMD_TPF,Psep))/Psep +1 +1; % +1 for background + 1 for remainder
end



for I = zero_MMD_TPF:Psep:4095
    
    [index_TPF,y] = find((fTHG_TPF_dot(:) >= I)&(fTHG_TPF_dot(:) <I+Psep)) ; % [x,y]為那些值的[row,column]的值,y已經被固定了，故不care。
    Sepa_TPF_dot(k) = struct('X_Axis' , index_TPF) ;
    k = k+1 ;
    
end


TPF_mean = 0;
fTHG_THG_mean = 0;
bTHG_THG_mean = 0;
fTHG_ST = 0;
bTHG_ST = 0;


for r = 1:Nsep
    
    TPF_mean(r) = mean(fTHG_TPF_dot(getfield(Sepa_TPF_dot(1,r) , 'X_Axis'))) ; % 求各分群中其螢光強度的平均值
    
end

for t = 1:Nsep
    
    fTHG_THG_mean(t) = mean (fTHG_THG_dot(getfield(Sepa_TPF_dot(1,t) , 'X_Axis'))) ; % 求相對各分群中螢光強度的有效三倍頻強度(平均值)
    bTHG_THG_mean(t) = mean (bTHG_THG_dot(getfield(Sepa_TPF_dot(1,t) , 'X_Axis'))) ; % 求相對各分群中螢光強度的有效三倍頻強度(平均值)
    fTHG_ST(t) = std (fTHG_THG_dot(getfield(Sepa_TPF_dot(1,t) , 'X_Axis')));
    bTHG_ST(t) = std (bTHG_THG_dot(getfield(Sepa_TPF_dot(1,t) , 'X_Axis')));
    
end




%% MMD Calibration 

  MMD_TPEF_eq = 'MMD = (TPEF - 390.757)/68.0094';
  MMD_mean = (TPF_mean - 390.757)/68.0094  ;%modified 20131028
        


MMD_mean(1) = 0; %set MMD = 0 for background data

%% Normalization (relative to background THG)
 
BG_fTHG = fTHG_THG_mean(1); 
BG_bTHG = bTHG_THG_mean(1);

fTHG_THGratio = fTHG_THG_mean / BG_fTHG; 
bTHG_THGratio = bTHG_THG_mean / BG_bTHG; 
bTHG_ST_ratio = bTHG_ST/ BG_bTHG; 



%% Regression (Phase 1)


%Set XY parameters
size_MMD_mean = size(MMD_mean);
size_MMD_mean = size_MMD_mean(2);

%Set data range for fitting


FitOne_MMD_mean = MMD_mean(1:size_FitOne);
FitOne_bTHG_THGratio = bTHG_THGratio(1:size_FitOne); 


%Choose positive THG data after bg correction

%FitOne_X_data : select (THGR - bgTHGR) > 0 
%FitOne_Y_data : THGR - bgTHGR

Count_FitOne = 0;
for i = 1: size_FitOne
    
    if FitOne_bTHG_THGratio(i) - bg_THGratio > 0
        Count_FitOne = Count_FitOne + 1;
      
        FitOne_X_data(Count_FitOne) = FitOne_MMD_mean(i);
        FitOne_Y_data(Count_FitOne) = FitOne_bTHG_THGratio(i) - bg_THGratio;
             
    end
end


FitOne_C = polyfit(log10(FitOne_X_data),log10(FitOne_Y_data),1);



%% Regression (Phase 2)

%Set data range for fitting
%Choose range of data

FitTwo_MMD_mean = MMD_mean(size_FitOne + 1: size_MMD_mean);
FitTwo_bTHG_THGratio = bTHG_THGratio(size_FitOne + 1: size_MMD_mean); 

%manipulate data for PL fitting
%Choose positive THG data after bg correction

%FitTwo_X_data : select (THGR - bgTHGR) > 0 
%FitTwo_Y_data : THGR - bgTHGR


FitTwo_X_data = MMD_mean(size_FitOne + 1: size_MMD_mean);
FitTwo_Y_data = bTHG_THGratio(size_FitOne + 1: size_MMD_mean ) - bg_THGratio;

%bg_THGratio = 1.06;

% POWER-LAW DISTRIBUTIONS fitting

FitTwo_C = polyfit(log10(FitTwo_X_data),log10(FitTwo_Y_data),1);






%% Intersection

x_intersection = 10.^((FitOne_C(2)-FitTwo_C(2))/(FitTwo_C(1)-FitOne_C(1)));
y_intersection = 10.^(FitOne_C(2))*x_intersection.^(FitOne_C(1)) + bg_THGratio;


%% ----- present -----

%% Plot normalized mean data  (THGR vs MMD)

%% Plot Figure 4.9 (a): mean+error (THGR vs MMD)

%% Plot Figure 4.9 (b):  x: 1-16 (loglog)


%Set parameters 

FitOne_X = linspace(min(FitOne_X_data),max(FitOne_X_data),100);
FitOne_Y = 10.^(FitOne_C(2))*FitOne_X.^(FitOne_C(1)); % not plus bgTHGR back yet. 


%Plot
figure,
loglog((FitOne_X_data),(FitOne_Y_data),'bo',(FitOne_X),(FitOne_Y),'--r') 


%Set title 
% Title = 'Empirical Equation';
% size_Title = 15 %set
% T = title(Title,'Fontsize',size_Title);

%Set legend 
size_legend = 19; %set

legend_Data = 'Data';
legend_FitOne_1_1 = '[THGratio - ';
legend_FitOne_1_2 = num2str(bg_THGratio);
legend_B_1_3 = '] = ';
legend_FitOne_1 = [legend_FitOne_1_1 legend_FitOne_1_2 legend_B_1_3];

legend_FitOne_2 = num2str(10.^(FitOne_C(2)));
legend_FitOne_3 = ['[MMD]','^{',num2str(FitOne_C(1)),'}']; 
legend_FitOne = [legend_FitOne_1 legend_FitOne_2 legend_FitOne_3];


L = legend(legend_Data,legend_FitOne);
set(L,'FontSize',size_legend);%set


%Set ticks
size_ticks = 20; 
set(gca,'FontSize',size_ticks) %Set font size to tick lables, the same as legend 

%Set label
size_label = 24; 
xlabel('Melanin Mass Density (mg/ml)','Fontsize',size_label)
ylabel('THG enhanced ratio','Fontsize',size_label)

%box
box on

%grid
grid on 

%Set axis
axis([1.5 12 0 8]) %20190308

 %% Plot Figure 4.9 (c):   x: 1-16 

%Set XY parameters

FitOne_plot_X = linspace(min(FitOne_MMD_mean),max(FitOne_MMD_mean),100);

FitOne_plot_Y = 10.^(FitOne_C(2))*FitOne_plot_X.^(FitOne_C(1)) + bg_THGratio;


%Hold 
hold on 

legend_Data = 'Data';
legend_FitOne_1_1 = '[THGratio ';
legend_FitOne_1_2 = '] = ';
legend_FitOne_1 = [legend_FitOne_1_1 legend_FitOne_1_2];

legend_FitOne_2 = num2str(10.^(FitOne_C(2)));
legend_FitOne_3 = ['[MMD]','^{',num2str(FitOne_C(1)),'}'];
legend_FitOne_4 = ' + ';
legend_FitOne_5 = num2str(bg_THGratio);

legend_FitOne = [legend_FitOne_1 legend_FitOne_2 legend_FitOne_3 legend_FitOne_4 legend_FitOne_5];

%%  Plot Figure 4.9 (d):   x: 17-75 (loglog)

%% Plot Figure 4.9 (e):   x: 17-75 
FitTwo_plot_X = linspace(min(FitTwo_MMD_mean),max(FitTwo_MMD_mean),100);
FitTwo_plot_Y = 10.^(FitTwo_C(2))*FitTwo_plot_X.^(FitTwo_C(1)) + bg_THGratio;


%Hold 
hold on 

legend_Data = 'Data';
legend_FitTwo_1_1 = '[THGratio ';
legend_FitTwo_1_2 = '] = ';
legend_FitTwo_1 = [legend_FitTwo_1_1 legend_FitTwo_1_2];

legend_FitTwo_2 = num2str(10.^(FitTwo_C(2)));
legend_FitTwo_3 = ['[MMD]','^{',num2str(FitTwo_C(1)),'}'];
legend_FitTwo_4 = ' + ';
legend_FitTwo_5 = num2str(bg_THGratio);

legend_FitTwo = [legend_FitTwo_1 legend_FitTwo_2 legend_FitTwo_3 legend_FitTwo_4 legend_FitTwo_5];

%% Plot Figure 4.9 (f):   whole range

% close all

%Plot
figure,
hold all
plot(MMD_mean,bTHG_THGratio,'ko') 
plot(FitOne_plot_X,FitOne_plot_Y,'--r','LineWidth',1.2);
plot(FitTwo_plot_X,FitTwo_plot_Y,'--g','LineWidth',1.2);


%Set title 
% Title = 'Empirical Equation';
% size_Title = 15 %set
% T = title(Title,'Fontsize',size_Title);

%Set legend 
size_legend = 19; %set
L = legend(legend_Data,legend_FitOne,legend_FitTwo);
set(L,'FontSize',size_legend);%set


%Set ticks
size_ticks = 20; 
set(gca,'FontSize',size_ticks) %Set font size to tick lables, the same as legend 

%Set label
size_label = 24; 
xlabel('Melanin Mass Density (mg/ml)','Fontsize',size_label)
ylabel('THG enhanced ratio','Fontsize',size_label)

%box
box on

%grid
grid on 

axis([0 60 0 28]) %20190308

%% Plot Figure 4.9 (g):   intersection point 

%% Plot Figure 4.9 (h):   empirical eq. 
