% This function takes input parameters 'Path_Network_Genes','Paired_File','Path_CNV' and 'Results_File'
% Returns the logical array of StructuralVariants for Tumor samples [Rows represent sorted list of Network Genes and Columns represents the order of CID in paired sample file]

function [StructuralVariants_LA]=Filter_Structural_Variants(Path_Network_Genes,Sheet_Filtered,Path_StructuralVariants,Results_File)
tic;

% Intersect network genes with CNVGenesList to find the index
% Network_Genes_List=Path_Network_Genes;
% NetworkGenes=readtable(Network_Genes_List);
% NetGenesCount=sum(NetworkGenes{:,2}~= "");

% Read Network Genes list
NetworkGenes=Path_Network_Genes;
T2 = readtable(NetworkGenes);
GenesList=table2cell(T2(:,1));
NGenesSum=sum(T2{:,2}~="");
NGenes=table2cell(T2(1:NGenesSum,2));

[SortedNGenes,GenesInd]=intersect(GenesList,NGenes);
Net_GNames=SortedNGenes;

% Get the paired samples from the sample sheet
% Note: If we make a separate function to get the logical array of CNV for
% our network genes and paired patients then we have to call functions
Filtered_Sample_Sheet=Sheet_Filtered;
% Number of columns in sheet differs in paired,4 unpaired_T,2 and
% OnlyTumor,3
if(size(Sheet_Filtered,2)==4)
    Tumor_CaseID=Filtered_Sample_Sheet(2:end,2);
elseif(size(Sheet_Filtered,2)==2)
    Tumor_CaseID=Filtered_Sample_Sheet(2:end,1);
else
    Tumor_CaseID=Filtered_Sample_Sheet(2:end,1);
end

% Get Mutations Data ColA and ColQ, mutated genes and SampleIDs stored in  'ColA1' and 'ColQ1'
StructuralVariants=Path_StructuralVariants;

if(strlength(StructuralVariants)==0)
    % Return empty logical array for SV of the same length
    StructuralVariants_LA=zeros(length(Net_GNames),length(Tumor_CaseID));
    display("Structural variants data is missing");
else
    
    SV_T=readtable(StructuralVariants);
    fieldnames(SV_T);
    ColA1 = SV_T.('Hugo_Symbol');
    ColD1=SV_T.('Tumor_Sample_Barcode');
    NewTable=[ColA1,ColD1];
    
    % Sorting the data so that it will help if we use intersect function
    SortedT=sortrows(NewTable,2);
    SV_GenesAll=SortedT(:,1);
    
    StructuralVariants_LA=zeros(length(Net_GNames),length(Tumor_CaseID));
    for i=1:length(Tumor_CaseID)
        indexfind=strfind(SortedT(:,2),char(Tumor_CaseID(i)));
        IndList=find(~cellfun(@isempty,indexfind));
        %All mutated genes per CID
        SV_Genes_CID=SV_GenesAll(IndList);
        [NetGenesSortedMutList,NetGenesInd]=intersect(SV_Genes_CID,Net_GNames);
        StructuralVariants_LA(:,i)=ismember(sort(Net_GNames),NetGenesSortedMutList);
    end
        % Add Headers in the Somatic Mutations logical array
    Tumor_Structural_Variants_LA=[['NetworkGenes';sort(Net_GNames)],[transpose(Tumor_CaseID);num2cell(StructuralVariants_LA)]];
    
    %% Write results data file %%
    Results_File_Path=cellstr(Results_File(1));
    UUID=cellstr(Results_File(2));
    
    % Tumor Raw Expression
    ResultsFilePath_T=strcat(char(Results_File_Path),'Tumor_Structural_Variants_Logical_Array_',char(UUID),'.xlsx');
    xlswrite(ResultsFilePath_T,Tumor_Structural_Variants_LA);
       
end
toc;
end