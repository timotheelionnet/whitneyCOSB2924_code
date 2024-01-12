%% data paths

% lambert et al Cell 2018 supplementary table listing the human TFs
% data in sheet 2
% first column is ensembl gene ID, second is gene name, third the DBD type
% and fourth a Yes/No TF qualifier (yes if thought to be sequence-specific TF, no for coactivators)
lambertPath = '/Users/lionnt01/Dropbox/articles/2023_COSB_chromatin_review/data/lambertCell2018/mmc2.xlsx';

% folder containing processed RNA-seq data for the deeply profiled Encode cell lines
% each dat file is in a separate tsv file.
% see the metadata.tsv file for the list of cell lines.
encodePath = '/Users/lionnt01/Dropbox/articles/2023_COSB_chromatin_review/data/encode';

%% load list of human TFs from Lambert supp table

lambertRaw = readtable(lambertPath,'Sheet',2);
lambertRaw = renamevars(lambertRaw, {'Var4'},{'isTF'});

% tflist is the list of ensembl IDs for all the SEQUENCE SPECIFIC human TFs
tfList = lambertRaw.ID(ismember(lambertRaw.isTF,{'Yes'}));

%% load encode RNA-seq data

tsvFilesRaw = dir(fullfile(encodePath, '*.tsv'));

% separate the data and metadatafiles
metaDataFile = tsvFilesRaw(ismember({tsvFilesRaw(:).name}','metadata.tsv'));
tsvFiles = tsvFilesRaw(~ismember({tsvFilesRaw(:).name}','metadata.tsv'));

% load the data
rnaSeqRaw = cell(numel(tsvFiles),1);
for i=1:numel(tsvFiles)
    rnaSeqRaw{i} = readtable(fullfile(tsvFiles(i).folder,tsvFiles(i).name),...
        'FileType','text','Delimiter','\t');
end

%% extract Ensembl gene IDs from the encode target_id column
% there are different formats for the different RNA-seq datasets (ugh)
    % first format:
    % the target ID is something like ''ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|processed_transcript|' in just the 'ENSG00000223972.5'
    % we just want the ENSG00000223972 (without the dot part) in order to match
    % the ID formats of the Lambert table
    % the tpm column is lower case.

    % second format (ignored for now): 
    % the ID is an ID like 10904
    % the tpm column is upper case
pattern = 'ENSG\d+';
extractENSG = @(x) regexp(x, pattern, 'match', 'once');
for i=1:numel(rnaSeqRaw)
    % deal with entries that follow format 1
    if ismember('tpm',rnaSeqRaw{i}.Properties.VariableNames)
        rnaSeqRaw{i}.ENSG = cellfun(extractENSG, rnaSeqRaw{i}.target_id, 'UniformOutput', false);
        rnaSeqRaw{i}.ENSG = cellfun(@(x) convertCharsToStrings(x), rnaSeqRaw{i}.ENSG, 'UniformOutput', false);
        rnaSeqRaw{i}.ENSG = string(rnaSeqRaw{i}.ENSG);
    end

end

% note: since the encode table lists all isoforms, there will be redundant ENSG
% entries.

%% 
% compute max tpms
maxTpm = 0;
for i=1:numel(rnaSeqRaw)
    if max(rnaSeqRaw{i}.tpm)>maxTpm
        maxTpm = max(rnaSeqRaw{i}.tpm);
    end
end

%% compute histogram of log(tpm+1) for all genes and for tfs only.
logTpmAllGenes = cell(numel(tsvFiles),1);
logTpmTFs = cell(numel(tsvFiles),1);

logTpmBinSize = 0.1;
histBins = 0:logTpmBinSize:log(maxTpm+1);

histLogTpmAllGenes = zeros(numel(rnaSeqRaw),numel(histBins));
histLogTpmTFs = zeros(numel(rnaSeqRaw),numel(histBins));

tfCounts = [];
allGeneCounts = [];
thresh1 = 0.8;
thresh2 = log(2)*2;
for i=1:numel(rnaSeqRaw)
    % sum the tpms for all entries that belong to the same gene (i.e. share
    % the same ENSG ID)
    % t is a table with columns ENSG | GroupCount | sum_tpm
    if ismember('ENSG',rnaSeqRaw{i}.Properties.VariableNames) ...
            && ismember('tpm',rnaSeqRaw{i}.Properties.VariableNames)
        t = groupsummary(rnaSeqRaw{i}, 'ENSG', 'sum', 'tpm');
    
        logTpmAllGenes{i} = log(t.sum_tpm + 1);
        logTpmTFs{i} = log(t.sum_tpm (ismember(t.ENSG(:)',tfList)) + 1);
        
        histLogTpmAllGenes(i,:) = hist(logTpmAllGenes{i}, histBins );
        histLogTpmTFs(i,:) = hist(logTpmTFs{i}, histBins );

        allGeneCounts = [allGeneCounts; [i, ...
            sum(logTpmAllGenes{i} >thresh1) ,...
            sum(logTpmAllGenes{i} >thresh2), numel(logTpmAllGenes{i})]];

        tfCounts = [tfCounts; [i, ...
            sum(logTpmTFs{i} >thresh1) ,...
            sum(logTpmTFs{i} >thresh2),numel(logTpmTFs{i}) ]];
    end
end

tfCounts = array2table(tfCounts,'VariableNames',{'datasetIdx',...
    'nGenesAboveThresh1','nGenesAboveThresh2','nGenesTotal'});

allGeneCounts = array2table(allGeneCounts,'VariableNames',{'datasetIdx',...
    'nGenesAboveThresh1','nGenesAboveThresh2','nGenesTotal'});

%% colect cell line names for the various indices
% load metadata file which has the key between file Accession name and cell
% line name
metaData = readtable(fullfile(metaDataFile.folder,metaDataFile.name),...
        'FileType','text','Delimiter','\t');

% connect each experiment with the cell line
for i=1:size(allGeneCounts)
    % find file name the data was taken
    [~,curFileName,~] = fileparts(tsvFiles(i).name)

    tfCounts.accession{i} = curFileName;
    allGeneCounts.accession{i} = curFileName;

    tfCounts.cellLine(i) = metaData.BiosampleTermName(...
        ismember(metaData.FileAccession,curFileName));
    allGeneCounts.cellLine(i) = metaData.BiosampleTermName(...
        ismember(metaData.FileAccession,curFileName));
end

% average the duplicate data for each cell line
tfCountsAvg = groupsummary(tfCounts,'cellLine',{'mean','std'},...
    {'nGenesAboveThresh1','nGenesAboveThresh2','nGenesTotal'});
allGeneCountsAvg = groupsummary(allGeneCounts,'cellLine',{'mean','std'},...
    {'nGenesAboveThresh1','nGenesAboveThresh2','nGenesTotal'});

%% plot histograms of tpm for all genes
figure('Name','histograms of RNA-seq expression');
hold on;
for i=1:numel(rnaSeqRaw)
    if sum(histLogTpmAllGenes(i,:)) ~=0
        plot(histBins,histLogTpmAllGenes(i,:)/sum(histLogTpmAllGenes(i,:)),...
            '-','color',[0.8,0.8,0.8],'DisplayName',['All genes, ',tfCounts.cellLine{tfCounts.datasetIdx==i}]);
    end
end

% plot histograms of tpm for tfs
for i=1:numel(rnaSeqRaw)
    if sum(histLogTpmTFs(i,:)) ~=0
        plot(histBins,histLogTpmTFs(i,:)/sum(histLogTpmTFs(i,:)),...
            '-','color',[0.6,0,0.8],'DisplayName',['TFs only, ',tfCounts.cellLine{tfCounts.datasetIdx==i}]);
    end
end

xlabel('RNAseq log(tpm+1)');
ylabel('Probability');
legend show;

%% plot counts of TF expressed per cell line
figure('Name','number of TFs expressed per cell line');
hold on;

tfCountsAvgSorted = sortrows(tfCountsAvg,'mean_nGenesAboveThresh2');
plot(tfCountsAvgSorted.mean_nGenesAboveThresh1,1:size(tfCountsAvgSorted,1),'o',...
    'MarkerEdgeColor',[0.8,0.8,0.8],'MarkerFaceColor',[0.8,0.8,0.8],...
    'DisplayName',['TFs only, inclusive expression threshold (log(tpm+1)>',...
    num2str(thresh1),')']);

plot(tfCountsAvgSorted.mean_nGenesAboveThresh2,1:size(tfCountsAvgSorted,1),'o',...
    'MarkerEdgeColor',[0.8,0.8,0.8],'MarkerFaceColor',[0.4,0.4,0.4],...
    'DisplayName',['TFs only, stringent expression threshold (log(tpm+1)>',...
    num2str(thresh2),')']);

xlabel('Number of TFs above expression threshold');
yticks(1:size(tfCountsAvgSorted,1));
yticklabels(tfCountsAvgSorted.cellLine);
legend show;

disp(['The median number of TFs expressed at log(tpm+1) > ',num2str(thresh2),...
    ' is ' ,num2str(median(tfCountsAvgSorted.mean_nGenesAboveThresh2))]);