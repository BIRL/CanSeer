%% This function takes parameters 'NormalRawExpressionFile','TumorRawExpressionFile' and 'Results_File'  
%% Returns the logical array of outliers for both Normal and Tumor samples [Rows represent sorted list of Network Genes and Columns represents the order of CID in paired sample file]
%% Returns detailed IQR computation file for both Normal and Tumor samples

function [N_Outliers_Ind,N_Outliers,Normal_IQR,T_Outliers_Ind,T_Outliers,Tumor_IQR]=Outliers_Detection_Using_Inter_Quartile_Range_Paired_Unpaired(NormalRawExpressionFile,TumorRawExpressionFile,Results_File)

tic;
% Get Normal Raw expression as the input parameter
NormalRawExpression=NormalRawExpressionFile;
[Rows,Cols]=size(NormalRawExpression);
NRawExp=NormalRawExpression(2:end,2:end);
NUM=cell2mat(NRawExp(:,:));
NCID=NormalRawExpression(1,2:end);
[NROWS,NCOLS]=size(NUM);
N_Outliers=zeros(NROWS,NCOLS);

% Variable Initialization for Normal IQR computation
N_Q1=[];
N_Q3=[];
N_IQR=[];
N_LowerBound=[];
N_UpperBound=[];
N_OutliersSum=[];
Normal=struct;
N_=struct;
N_Count_Outliers=0;

% Get Network Gene Names from the expression data file 
Net_Genes=NormalRawExpression(:,1);

% Get Tumor Raw expression as the input parameter
TumorRawExpression=TumorRawExpressionFile;
TRawExp=TumorRawExpression(2:end,2:end);
NUM1=cell2mat(TRawExp(:,:));
CCID=TumorRawExpression(1,2:end);
[TROWS,TCOLS]=size(NUM1);
T_Outliers=zeros(TROWS,TCOLS);

% Variable Initialization for Tumor IQR computation
T_Q1=[];
T_Q3=[];
T_IQR=[];
T_LowerBound=[];
T_UpperBound=[];
T_OutliersSum=[];
Cancer=struct;
C_=struct;
T_Count_Outliers=0;

% Compute this loop for number of network genes 'NROWS'
for i=1:NROWS
    % Get expression values from all the from each gene row 'i' and all patients columns ':'
    Out_temp= NUM(i,:);
    temp=sort(Out_temp);
    
    % Compute Q1, Q3 and IQR of every network gene across all paired normal samples
    NN_Q1= median(temp(find(temp<median(temp)))); %Compute 25th percentile
    NN_Q3= median(temp(find(temp>median(temp)))); %Compute 75th percentile
    NN_IQR=NN_Q3-NN_Q1;
    
    % Use 'Out_temp'unsorted patients expression data to compute correct indexes of outliers
    NNTemp= Out_temp < NN_Q1 - (1.5 * NN_IQR);
    NNTemp1= Out_temp > NN_Q3 + (1.5 * NN_IQR);
    xx= NNTemp + NNTemp1;
    N_Outliers(i,:)= xx;
    Nsum=length(find(xx));
    
    % IQR Values stored in cell arrays
    N_Q1= [N_Q1,num2cell(NN_Q1)];
    N_Q3= [N_Q3,num2cell(NN_Q3)];
    N_IQR= [N_IQR,num2cell(NN_IQR)];
    N_LowerBound=[N_LowerBound,num2cell(NN_Q1 - (1.5 * NN_IQR))];
    N_UpperBound=[N_UpperBound,num2cell(NN_Q3 + (1.5 * NN_IQR))];
    N_OutliersSum=[N_OutliersSum,num2cell(Nsum)];
    
    % Get expression values from all the from each gene row 'i' and all patients columns ':'
    Out_temp1=NUM1(i,:);
    temp1=sort(Out_temp1);
    % Compute Q1, Q3 and IQR to find lower and upper bound to filter
    % outliers of every  gene across normal paired samples and then
    % cancer samples separately
    TT_Q1= median(temp1(find(temp1<median(temp1))));
    TT_Q3= median(temp1(find(temp1>median(temp1))));
    TT_IQR= TT_Q3 - TT_Q1;
    TTemp= Out_temp1 < TT_Q1 - (1.5 * TT_IQR);
    TTemp1= Out_temp1 > TT_Q3 + (1.5 * TT_IQR);
    yy= TTemp + TTemp1;
    T_Outliers(i,:)= yy;
    Tsum=length(find(yy));
    
    % IQR Values stored in cell arrays
    T_Q1=[T_Q1,num2cell(TT_Q1)];
    T_Q3=[T_Q3,num2cell(TT_Q3)];
    T_IQR=[T_IQR,num2cell(TT_IQR)];
    T_LowerBound=[T_LowerBound,num2cell(TT_Q1 - (1.5 * TT_IQR))];
    T_UpperBound=[T_UpperBound,num2cell(TT_Q3 + (1.5 * TT_IQR))];
    T_OutliersSum=[T_OutliersSum,num2cell(Tsum)];
    
    % Find index of outliers of gene across Normal Samples
    N_tempo=xx;
    Ind=find(N_tempo);
    Normal(i).GeneInd=num2cell(Ind);
    N_(i).Gene=Ind;
    [R,C]=size(N_(i).Gene);
    N_tempo=C;
    if(N_tempo>N_Count_Outliers)
        N_Count_Outliers=N_tempo;
    end
    % Find index of outliers of gene across Cancer Samples
    T_tempo=yy;
    Ind=find(T_tempo);
    Normal(i).GeneInd=num2cell(Ind);
    C_(i).Gene=Ind;
    [R,C]=size(C_(i).Gene);
    T_tempo=C;
    if(T_tempo>T_Count_Outliers)
        T_Count_Outliers=T_tempo;
    end   
end
% Concatenate IQR values with Expression file
Header={'Quartile1', 'Quartile3', 'IQR', 'LowerBound', 'UpperBound','SumOfOutliers'};
Normal=[N_Q1;N_Q3; N_IQR; N_LowerBound; N_UpperBound;N_OutliersSum];
Normal_IQR=horzcat(NormalRawExpressionFile,(vertcat(Header,transpose(Normal))));

Tumor=[T_Q1;T_Q3; T_IQR;T_LowerBound; T_UpperBound;T_OutliersSum];
Tumor_IQR=horzcat(TumorRawExpressionFile,(vertcat(Header,transpose(Tumor))));

% Index of the outliers stored in matrix and '0' index in the matrix is not
% associated with any patient so ignore '0' values in the matrix
%Ignore zeros in this table
N_Outliers_Ind=zeros(NROWS,N_Count_Outliers);
T_Outliers_Ind=zeros(NROWS,T_Count_Outliers);
for i=1: NROWS
    [R,C]=size(N_(i).Gene);
    N_Outliers_Ind(i,1:C)=N_(i).Gene;
    [R1,C1]=size(C_(i).Gene);
    T_Outliers_Ind(i,1:C1)=C_(i).Gene;
end

% Include Headers in the Outliers logical array
Header_Normal= [NCID;num2cell(N_Outliers)];
Header_Tumor=[CCID;num2cell(T_Outliers)];
Normal_Outliers_LA=[Net_Genes,Header_Normal];
Tumor_Outliers_LA=[Net_Genes,Header_Tumor];

%% Write results data file %%
Results_File_Path=cellstr(Results_File(1));
UUID=cellstr(Results_File(2));
% Normal IQR and logical array
ResultsFilePath_N=strcat(char(Results_File_Path),'Normal_Outliers_Logical_Array_IQR_',char(UUID),'.xlsx');
xlswrite(ResultsFilePath_N,Normal_Outliers_LA);
ResultsFilePath_NN=strcat(char(Results_File_Path),'Normal_IQR_',char(UUID),'.xlsx');
xlswrite(ResultsFilePath_NN,Normal_IQR);
% Tumor IQR and logical array
ResultsFilePath_T=strcat(char(Results_File_Path),'Tumor_Outliers_Logical_Array_IQR_',char(UUID),'.xlsx');
xlswrite(ResultsFilePath_T,Tumor_Outliers_LA);
ResultsFilePath_TT=strcat(char(Results_File_Path),'Tumor_IQR_',char(UUID),'.xlsx');
xlswrite(ResultsFilePath_TT,Tumor_IQR);

toc;
end