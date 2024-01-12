
% data from "Integrated spatial genomics reveals global architecture of
% single nuclei" Takei et al, Nature 2021
% positions downloaded from zenodo at https://zenodo.org/records/3735329
% locus Key is supplementary table 1 at https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-03126-2/MediaObjects/41586_2020_3126_MOESM2_ESM.xlsx

% i only use the 25 kbp resolution dataset. 
% the dataset contains the locations of 60 adjacent regions on each of the chromosomes.
% so the regions are chr1-1 chr1-2, chr1-3, ..., chr1-60, chr2-1, chr2-2,
% ... chr2-60 etc.

posRep1Path = '/Users/lionnt01/Dropbox/articles/2023_COSB_chromatin_review/data/DNAseqFISH+ 2/DNAseqFISH+25kbloci-E14-replicate1.csv';
posRep2Path = '/Users/lionnt01/Dropbox/articles/2023_COSB_chromatin_review/data/DNAseqFISH+ 2/DNAseqFISH+25kbloci-E14-replicate2.csv';
locusKeyPath = '/Users/lionnt01/Dropbox/articles/2023_COSB_chromatin_review/data/takeiNature2021_SuppTable1_41586_2020_3126_MOESM2_ESM.xlsx';

xPixToNm = 103;
yPixToNm = 103;
zPixToNm = 250;
addpath('cbrewer');
%% load locus key
% format is:
% Region ID, a unique identifier for each locus
% Name: a name following the format 'chr1-HighRes-#31' where "31" in thisn example a number 1-60 that encodes the location of the spot on a chromosome)
% Channel: fluorescence channel is always 3 for this experiment
% Chrom: chromosome ID in the format chr1, chr2, etc
% Start: start coordinate in bp
% End: End coordinate in bp
% ChromID: chromosome ID as a number
locKeyRaw = readtable(locusKeyPath);

%% load data replicates
% format is:
% fov: field of view ID
% cell: cell ID
% region ID: the ID of the region (number 1-60; needs to cbe combined with
% the chromID to obtain the locus unambiguously)
% x y z coordinates in pixel units (pixel sizes are 103 nm for x and y, and 250 nm for z)
% dot intensity: (?) intensity of the spot in the 1-60 hybridizations where
    % each chromosome is walked in parallel
% chr1-X intensity: intensity of the chromosome painting channel used to
% assign chromosome identity.
% chromID: chromosome ID assignment based on the chr intensities.
% labelID: allele ID -- can be -1 (I assume problem alleles), 0, 1, 2, 3
pos1raw = readtable(posRep1Path);
pos2raw = readtable(posRep2Path);

%% reformat locusKey
locKey = removevars(locKeyRaw,{'Channel','RegionID'});
locKey.Pos = (locKey.Start+locKey.End)/2;
locKey = removevars(locKey,{'Start','End','Chrom'});

% Extract the numeric part from each string using cellfun and regexp
numericParts = cellfun(@(x) regexp(x, '(?<=#)\d+', 'match'), locKey.('Name'), 'UniformOutput', false);

% Since regexp returns a cell array, you might get nested cells. Flatten them.
numericParts = cellfun(@(x) x{1}, numericParts, 'UniformOutput', false);

% Update the table column
locKey.('Name') = numericParts;

% convert the column to numeric
locKey.('Name') = str2double(locKey.('Name'));

% rename to regionID
locKey = renamevars(locKey,'Name','regionID');

%% reformat distance data
pos1 = reformatDistTable(pos1raw,xPixToNm,yPixToNm,zPixToNm);
pos2 = reformatDistTable(pos2raw,xPixToNm,yPixToNm,zPixToNm);

%% plot the distribution of genomic distances profiled
chrList = unique(locKey.ChromID);
genomicDistFull = [];
for i=1:numel(chrList)
    tmp = locKey.Pos(locKey.ChromID == chrList(i));
    dtmp = abs(repmat(tmp,1,60) - repmat(tmp',60,1));
    genomicDistFull = [genomicDistFull;dtmp(:)];
end

genomicDistUnique = unique(genomicDistFull);
genomicDistUniqueN = zeros(size(genomicDistUnique));
for i=1:numel(genomicDistUnique)
    genomicDistUniqueN(i) = sum(genomicDistFull == genomicDistUnique(i));
end

figure('Name','Genomic Distances Probed');
plot(genomicDistUnique/1000,genomicDistUniqueN,'o','LineWidth',2);
xlabel('distance in kbp');
ylabel('number of loci pairs');

%% generate an example plots of the x,y coordinates of three chromosomes picked at random

% pick a replicate at random
repToPlot = randsample(1:2,1);
if repToPlot == 1
    p = pos1;
else
    p = pos2;
end

% pick a fov at random
fovList = unique(p.fov);
fovToPlot = randsample(fovList,1);
disp(['Found ',num2str(numel(fovList)),' fovs in replicate ',num2str(repToPlot),...
    '; plotting results from fov ',num2str(fovToPlot),', picked at random.']);
p2 = p(p.fov == fovToPlot,:);
cellList = unique(p2.cellID);

% pick a cell at random
cellToPlot = randsample(cellList,1);
p2 = p2(p2.cellID == cellToPlot,:);
disp(['Found ',num2str(numel(cellList)),' cells',...
    '; plotting results from cell ',num2str(cellToPlot),', picked at random.']);

%pick three chromosomes at random
chromToPlot = randsample(1:20,3);
p2 = p2(ismember(p2.chromID,chromToPlot),:);
disp(['plotting results from chromosomes ',num2str(chromToPlot),', picked at random.']);

% plot x,y coordinates
figure('Name',['x,y rep ',num2str(repToPlot),'; fov ',num2str(fovToPlot),...
    '; cell ',num2str(cellToPlot),'; chromosomes ',num2str(chromToPlot)]);
hold on;
for i=1:numel(chromToPlot)
    ptmp = p2(p2.chromID == chromToPlot(i),:);
    labelsToPlot = unique(ptmp.labelID);
    for j=1:numel(labelsToPlot)
        ptmp2 = ptmp(ptmp.labelID == labelsToPlot(j),:);
        plot(ptmp2.x,ptmp2.y,'linewidth', 2,'DisplayName',...
            ['Chrom ',num2str(chromToPlot(i)),'; label ',num2str(labelsToPlot(j))]);
    end
end
legend show;

%% prepare to compute all distances between pairs

% generate map where regionIDs are keys and chromosome regions coordinates 
% are values
regionKeys = cell(numel(chrList),1);
for i=1:numel(chrList)
    regionKeys{i} = containers.Map(...
        locKey.regionID(locKey.ChromID == chrList(i)),...
        locKey.Pos(locKey.ChromID == chrList(i)));
end

% compute the number of pairs in each allele to enable preallocating the size of the distance table
nr = 0;
varNames = {'rep', 'fov', 'cell', 'chr', 'allele','nPairs','cumsum'}; 
nPairsTable = table('Size', [0, length(varNames)],...
        'VariableTypes', repmat({'double'}, 1, length(varNames)), ...
        'VariableNames', varNames);  
for i1 = 1:2
    if i1 == 1
        p = pos1;
    else
        p = pos2;
    end
    fovList = unique(p.fov);
    for i2 = 1:numel(fovList)
        p2 = p(p.fov == fovList(i2),:);
        cellList = unique(p2.cellID);
        for i3 = 1:numel(cellList) 
            p3 = p2(p2.cellID == cellList(i3),:);
            for i4 = 1:numel(chrList)
                p4 = p3(p3.chromID == chrList(i4),:);
                alleleList = setdiff(unique(p4.labelID),-1);
                for i5 = 1:numel(alleleList)
                    p5 = p4(p4.labelID == alleleList(i5),:);
                    np = size(p5,1)*(size(p5,1)-1)/2;
                    curNp = table('Size', [1, length(varNames)],...
                        'VariableTypes', repmat({'double'}, 1, length(varNames)), ...
                        'VariableNames', varNames);  
                    curNp.rep(:) = i1;
                    curNp.fov(:) = fovList(i2);
                    curNp.cell(:) = cellList(i3);
                    curNp.chr(:) = chrList(i4);
                    curNp.allele(:) = alleleList(i5);
                    curNp.nPairs = np;
                    nPairsTable = [nPairsTable;curNp];
                end
            end
        end
    end
end
nPairsTable.cumsum = cumsum(nPairsTable.nPairs);
disp(['Preallocating distance table with ',num2str(sum(nPairsTable.nPairs)),' rows...']);
%% compute all distances between pairs
% create master table of distances
varNames = {'rep', 'fov', 'cell', 'chr', 'allele','regionID1','regionID2',...
    'genomicPos1','genomicPos2','genomicDist','physDist'}; 
% distTable = table('Size', [sum(nPairsTable.nPairs), length(varNames)],...
%         'VariableTypes', repmat({'double'}, 1, length(varNames)), ...
%         'VariableNames', varNames);  

x = zeros(sum(nPairsTable.nPairs),11);
distTable = table('Size', [1000000, length(varNames)],...
         'VariableTypes', repmat({'double'}, 1, length(varNames)), ...
         'VariableNames', varNames);  
ctr = 0;
for i1 = 1:2
%for i1 = 1:1
    if i1 == 1
        p = pos1;
    else
        p = pos2;
    end
    disp(['replicate ',num2str(i1),'...']);
    fovList = unique(p.fov);
    for i2 = 1:numel(fovList)
        % build a table for the fov
        fovPairs = nPairsTable(nPairsTable.rep == i1 ...
            & nPairsTable.fov == fovList(i2),:);
        nPairsInFov = sum(fovPairs.nPairs);
        fovStart = fovPairs.cumsum(1) - fovPairs.nPairs(1);
    
        fovT = table('Size', [nPairsInFov, length(varNames)],...
                'VariableTypes', repmat({'double'}, 1, length(varNames)), ...
                'VariableNames', varNames); 
        %for i2 = 1:1
        disp(['fov ',num2str(fovList(i2)),'/',num2str(numel(fovList)),'...']);
        p2 = p(p.fov == fovList(i2),:);
        cellList = unique(p2.cellID);
        disp(['Found ',num2str(numel(cellList)),' cells.']);
        
        for i3 = 1:numel(cellList)    
            % for i3 = 1:1
            tic;
            p3 = p2(p2.cellID == cellList(i3),:);
            for i4 = 1:numel(chrList)
            %for i4 = 1:1
            %disp(['chr ',num2str(i4),'/',num2str(numel(chrList)),'...']);
                p4 = p3(p3.chromID == chrList(i4),:);

                % allele IDs are 0,1,2,3 (ignoring the -1 values that I
                % assume correspond to problem alleles)
                alleleList = setdiff(unique(p4.labelID),-1);
                for i5 = 1:numel(alleleList)
                %disp(['allele ',num2str(alleleList(i5)),' ...']);
                    p5 = p4(p4.labelID == alleleList(i5),:);
                    n = size(p5,1);
                    % compute distance matrix
                    d = sqrt ( (repmat(p5.x,1,n) - repmat(p5.x',n,1)).^2 ...
                        + (repmat(p5.y,1,n) - repmat(p5.y',n,1)).^2 ...
                        + (repmat(p5.z,1,n) - repmat(p5.z',n,1)).^2 );

                    % generate similar matrices for regionIDs of i and j
                    ri = repmat(p5.regionID,1,n);
                    rj = repmat(p5.regionID',n,1);

                    % generate similar matrices for genomic coordinates of i and j
                    gi = cell2mat(values(regionKeys{i4},num2cell(ri)));
                    gj = cell2mat(values(regionKeys{i4},num2cell(rj)));

                    % generate similar matrix for genomic separation
                    gDist = gi-gj;
                    
                    % build a list of the linear indices of the bottom
                    % triangle of the distance matrix (excluding the
                    % diagonal) indices correspond to the order ([r2,c1],[r3,c1],[r4,c1],...,[r3,c2],[r4,c2],...)
                    % where the distance are ( r-c)
                    lowerTri = tril(d, -1);
                    l = find(lowerTri);
                    
                    % build a table for all the distances in the current allele
                    curT = table('Size', [max(0,((n-1)*n)/2), length(varNames)],...
                        'VariableTypes', repmat({'double'}, 1, length(varNames)), ...
                        'VariableNames', varNames);   
                    curT.rep(:) = i1;
                    curT.fov(:) = fovList(i2);
                    curT.cell(:) = cellList(i3);
                    curT.chr(:) = chrList(i4);
                    curT.allele(:) = alleleList(i5);
                    curT.genomicPos1(:) = gi(l);
                    curT.genomicPos2(:) = gj(l);
                    curT.regionID1(:) = ri(l);
                    curT.regionID2(:) = rj(l);
                    curT.genomicDist(:) = gDist(l);
                    curT.physDist(:) = d(l);
                    
                    idxAllele = fovPairs.cell == cellList(i3) ...
                        & fovPairs.chr == chrList(i4) ...
                        & fovPairs.allele == alleleList(i5);

                    idxStart = fovPairs.cumsum(idxAllele)- fovPairs.nPairs(idxAllele) - fovStart +1;
                    idxEnd = fovPairs.cumsum(idxAllele) - fovStart;
                    fovT(idxStart:idxEnd,:) = curT;
                    %ctr = ctr+size(curT,1);
                end
            end 
            t = toc;
            disp(['cell ',num2str(cellList(i3)),'/',num2str(numel(cellList)),...
                ' processed ',num2str(size(p3,1)),' spot positions in ',num2str(t),' seconds.']);            
        end
        distTable(fovStart+1:fovStart+nPairsInFov,:) = fovT;
    end
end
disp('Done.');

%% validation of the calculation with random picks
nRowsToCheck = 3;
idx = randsample(size(distTable,1),nRowsToCheck);
for i=1:numel(nRowsToCheck)
    disp('   ');
    disp('*************************************');
    disp(['Checking row ',num2str(idx(i)),'...'])
    r = distTable(idx(i),:)
    curRep = r.rep;
    curFov = r.fov;
    curCell = r.cell;
    curChrom = r.chr;
    curAllele = r.allele;
    disp(['   rep ',num2str(curRep),', fov ',num2str(curFov),...
        ', cell ',num2str(curCell),', chrom  ',num2str(curChrom),...
        ', allele ',num2str(curAllele)]);
    if curRep == 1
        p = pos1;
    else
        p = pos2;
    end
    g1 = r.genomicPos1;
    g2 = r.genomicPos2;
    regionID1 = locKey.regionID( locKey.ChromID == curChrom ...
        & locKey.Pos == g1);
    regionID2 = locKey.regionID( locKey.ChromID == curChrom ...
        & locKey.Pos == g2);

    p1 = p( p.fov == curFov & p.cellID == curCell & p.chromID == curChrom ...
        & p.labelID == curAllele & p.regionID == regionID1,:);

    p2 = p( p.fov == curFov & p.cellID == curCell & p.chromID == curChrom ...
        & p.labelID == curAllele & p.regionID == regionID2,:);
    n1 = size(p1,1);
    n2 = size(p2,1);
    disp(['    Found ',num2str(n1),' spot(s) for the locus''s first anchor and ',...
        num2str(size(p2,1)),' for the second anchor']);
    for i1=1:n1
        for i2 = 1:n2
        if (n1>1 || n2 > 1)
            disp(['    Anchors ',num2str(i1),' and ',num2str(i2),':']);
        end
        disp(['    Region 1 ID ',num2str(regionID1),' at ',num2str(g1),...
            '; x = ',num2str(p1.x(i1)),'; y = ',num2str(p1.y(i1)),'; z = ',num2str(p1.z(i1))]);
        disp(['    Region 2 ID ',num2str(regionID2),' at ',num2str(g2),...
            '; x = ',num2str(p2.x(i2)),'; y = ',num2str(p2.y(i2)),'; z = ',num2str(p2.z(i2))]);
        disp(['    Distance from distTable: ',num2str(r.physDist)]);
        d = sqrt( ( p1.x(i1) - p2.x(i2))^2 + ( p1.y(i1) - p2.y(i2))^2  + ( p1.z(i1) - p2.z(i2))^2  );
        disp(['    Validation distance calculation: ',num2str(d)]);
        disp(' ');
        end
    end
end

%% plot cfs of physical distances for all genomic separations probed.
cdfMin = min(distTable.physDist);
cdfMax = max(distTable.physDist);
distBinSize = 10;
distData = cell(numel(genomicDistUnique),1);
for i=1:numel(genomicDistUnique)
    distData{i} = distTable.physDist(distTable.genomicDist == genomicDistUnique(i));
    cdf{i} = cumsum(hist(distData{i},cdfMin:distBinSize:cdfMax))/numel(distData{i});
    medianDist(i) = median(distData{i});
    meanDist(i) = mean(distData{i});
end

figure('Name','cdfs');
hold on;
for i=2:numel(genomicDistUnique)
    plot(cdfMin:distBinSize:cdfMax,cdf{i},'Linewidth',2,...
        'DisplayName',['genomic Dist: ',num2str(genomicDistUnique(i))]);
end
xlabel('Physical Distance (nm)');
ylabel('CDF');

%% plot median physical distance as a function of genomic separation 
fMedian = figure('Name','median dist');
hold on;
plot(genomicDistUnique(2:end)/1000, medianDist(2:end),'-','Linewidth',8,'Color',[0.8,0.8,0.8]);
xlabel('Genomic separation (kbp)');
ylabel('Median physical Distance (nm)');
xlim([0,1000]);
fontsize(16,'points');
aMedian = gca;

%% plot mean physical distance as a function of genomic separation
figure('Name','mean dist');
plot(genomicDistUnique(2:end)/1000, meanDist(2:end),'-','Linewidth',8,'Color',[0.8,0.8,0.8]);
xlabel('Genomic separation (kbp)');
ylabel('Mean physical Distance (nm)');
xlim([0,1000]);
fontsize(16,'points');

%% display how many loci pairs were used
%compute all unique pairs chrID, regionID1, regionID2
x = unique(distTable(:,[4,6,7]),'rows');
disp(['Found ',num2str(size(x,1)),' loci pairs in full table (including doublets)']);
x = table2array(x);
x = [x(:,1),sort(x(:,2:3),2,'ascend')];
x = unique(x,'rows');    
disp(['Found ',num2str(size(x,1)),' loci pairs in full table after excluding doublets']);
%% add data from other sources
% formatted as 
% genomic separation in kbp | median distance in nm | transcription state
% (0 for inactive, 1 for active).

% Sox2 data from Alexander et al elife 2019 (LacO-cuO)
% https://elifesciences.org/articles/41769
sox2Alexander = [   117,255,0;... % separation from text
                    117,263,1]; % mean distance estimated by digitizing Fig. 6D - note that the earlier figures show larger separations. Different genotype backgrounds?

% Sox2 data from Platania et al Biorxiv 2023 (ANCH1-ANCH3)
% https://www.biorxiv.org/content/10.1101/2023.04.25.538222v1
sox2Platania = [    116,145,0;... % coordinates non obvious from text so estimated separation from genomic map in fig. 1A
                    116,166,1]; % data from fig. 2D comparing WT to the Sox2-SCR DeltaSRR107+111

% SHH from Benabdallah et al Mol Cell 2019
% https://pubmed.ncbi.nlm.nih.gov/31494034/
% SBE6 is 100 kbp away from SHH; SBE4 is 350; SBE2/3 450; ZRS 800
shhBenabdallah = [  100,327,0;...
                    100,423,1;... % separation from Fig. 1A, distances from 1B
                    350,342,0;...
                    350,464,1;...
                    450,408,0;...
                    450,480,1;...
                    800,398,0;...
                    800,408,1];

% synthetic loci in Drosophila from Bruckner et al Science 2023
% https://www.science.org/doi/10.1126/science.adf5568
flyBruckner = [     58,516,0; ... % separations from Fig. 2A, distances from figure 2B
                    58,377,1; ...
                    82,668,0;...
                    82,451,1;...
                    88,602,0;...
                    88,381,1;...
                    149,722,0;...
                    149,378,1;...
                    190,797,0;...
                    190,354,1;...
                    595,983,0;...
                    595,329,1];

% BX-C locus in drosophila from Mateo et al Nature 2019 
% https://www.nature.com/articles/s41586-019-1035-4
% distances iand separations are estimated from  extended data figure 7a,
% using the enhancer center as a reference point.
% abx UBX 53 kbp
% BRE UBX 35 kbp
% bxd UBX 12 kbp
% pbx UBX 35 kbp
% iab2 Abd-A 20 kbp
% iab3 Abd-A 15 kbp
% iab4 Abd-A 35 kbp
% iab5 Abd-B 59 kbp 
% iab6 Abd-B 44 kbp 
% iab7 Abd-B 22 kbp
% iab8 Abd-B 14 kbp

flyMateo = [53,359,0;...
            53,345,1;...
            35,352,0;...
            35,326,1;...
            12,326,0;...
            12,289,1;...
            35,343,0;...
            35,317,1;...
            20,308,0;...
            20,275,1;...
            15,343,0;...
            15,299,1;...
            35,337,0;...
            35,313,1;...
            59,358,0;...
            59,338,1;...
            44,322,0;...
            44,301,1;...
            22,313,0;...
            22,294,1;...
            14,285,0;...
            14,275,1];

% Ohishi Biorxiv 2023 https://www.biorxiv.org/content/10.1101/2023.11.27.568629v1
% figure 4c genomic separation based on labels,
% distance very approximate due to figure resolution
% including only the -45 and 60 kbp superenhancers.
nanogOhishi = [ 45,313,0;...
                45,313,1;...
                60,375,0;...
                60,375,1];

% Barinov et al Arxiv 2020. https://arxiv.org/abs/2012.15819
% even-skipped - note that the "inactive state" include the interstripe
% region where enhancers and the gene are inactive, plus stripes where the
% enhancer in question is inactive but the gene is expressing due to
% another enhancer.
% genomic separations from fig. 1A (using enhancer center to TSS). distances from fig. 2B
% eve-E5 5.8kbp
% eve-E4+6
% first inactive state for each pair is the interstripe
flyBarinov = [  5.8,165,0;...
                5.8,208,0;...
                5.8,186,0;...
                5.8,185,0;...
                5.8,202,0;...
                5.8,164,1;...
                3.1,165,0;...
                3.1,214,0;...
                3.1,199,0;...
                3.1,204,0;...
                3.1,182,1;...
                3.1,179,1];

% Chen et al Molecular Cell 2023. https://pubmed.ncbi.nlm.nih.gov/36996812/
% Sox9 separations E1.45: 1450 kbp E1.35: 1350 kbp, E1.25: 1250 kbp
% separations based on name of enhancers. Distances measured from Fig 1G
% (Mean distance)

sox9Chen = [1450,494,0;...
            1450,348,1;...
            1350,416,0;...
            1350,319,1;...
            1250,445,0;...
            1250,358,1];

% plot median physical distance as a function of genomic separation 
fMedian = figure('Name','median dist');
hold on;
plot(genomicDistUnique(2:end)/1000, medianDist(2:end),'-','Linewidth',8,'Color',[0.8,0.8,0.8]);
xlabel('Genomic separation (kbp)');
ylabel('Median physical Distance (nm)');
xlim([0,1500]);
fontsize(16,'points');
aMedian = gca;

dataList = {sox2Alexander,sox2Platania,shhBenabdallah,flyBruckner,flyMateo,nanogOhishi,flyBarinov,sox9Chen};
dataLabels = {'sox2Alexander','sox2Platania','shhBenabdallah','flyBruckner','flyMateo','nanogOhishi','flyBarinov','sox9Chen'};

% plot links between on and off 
colorLink = [0.9,0.9,0.9];
for i=1:numel(dataList)
    distList = unique(dataList{i}(:,1)); % list unique separations
    for j=1:numel(distList)
        p = dataList{i}(dataList{i}(:,1)==distList(j),:);
        plot(aMedian,p(:,1),p(:,2),'-','color',colorLink,'Linewidth',2);
    end
end

% plot datapoints
colorOn = cbrewer('qual','Set1',numel(dataList)); % generating a larger palette to avoid the extreme colors
colorOff = [0.9,0.9,0.9];
for i=1:numel(dataList)
    for j=1:size(dataList{i},1)
        switch dataList{i}(j,3)
            case 0
                curColor = colorOff;
            case 1
                curColor = colorOn(i,:);
        end
        plot(aMedian,dataList{i}(j,1),dataList{i}(j,2),'o',...
            'MarkerFaceColor',curColor,'MarkerEdgeColor',curColor,...
            'MarkerSize',8,'DisplayName',dataLabels{i});
    end
end

%% Compute Cdfs for each of the pairs in the shortest distance bin (25 kbp)
% list all unique loci pairs separated by 25 kbp
pairIdx25 = unique(distTable(:,[4,6,7,10]),'rows');
pairIdx25 = pairIdx25(pairIdx25.genomicDist==25000,:);

% initialize cdf table
cdf25Bin = 10;
cdf25max = 3000;
x25 = 0:cdf25Bin:cdf25max;
cdf25Table = zeros(size(pairIdx25,1),numel(x25));
median25 = zeros(size(pairIdx25,1),1);
% loop through loci pairs separated by 25 kbp
for i=1:size(pairIdx25,1)

    % collect all measurements for a given pair of loci
    d = distTable( distTable.chr == pairIdx25.chr(i) ...
        & distTable.regionID1 == pairIdx25.regionID1(i) ...
        & distTable.regionID2 == pairIdx25.regionID2(i) ...
        ,:);
    
    % build cdf
    cdf25Table(i,:) = cumsum( hist(d.physDist(:),x25) )/size(d,1);
    median25(i,1) = median(d.physDist(:));
end

%% plot cdfs 25 kbp
% sort cdfs by their median distance
[median25sorted,idx25Sorted] = sort(median25,'ascend');
cdf25TableSorted = cdf25Table(idx25Sorted,:);

figure('Name','cdfs of loci separated by 25 kbp');
hold on;

percentBins = [1,2:2:98,99];
c25maps = cbrewer('div','PiYG',numel(percentBins));
for i=1:numel(percentBins)
    j = ceil( numel(median25sorted)*percentBins(i)/100);
    plot(x25,cdf25TableSorted(j,:),'Linewidth',2,'Color',c25maps(i,:),...
        'DisplayName',[num2str(percentBins(i)),'th percentile ']);
end

% overaly the cdf of the median distances (crowds the fig a bit)
% plot(median25sorted,(1:numel(median25sorted))/numel(median25sorted),...
%     'Linewidth',4,'Color',[0.6,0.6,0.6],...
%     'DisplayName','cdf of median distances for all pairs');

xlim([0,1000]);
xticks(0:200:1000);
xlabel('distance (nm)');
ylabel('cdf');
fontsize(16,'points');

colormap(c25maps);
colorbar;
%%
%****************************************************************************
% SUB FUNCTIONS
%****************************************************************************


% reformat the distance table
function d = reformatDistTable(draw,xPixToNm,yPixToNm,zPixToNm)
    d = draw;
    d = removevars(d,{'channel','dot_intensity','chr1_intensity','chr2_intensity',...
        'chr3_intensity','chr4_intensity','chr5_intensity','chr6_intensity',...
        'chr7_intensity','chr8_intensity','chr9_intensity','chr10_intensity',...
        'chr11_intensity','chr12_intensity','chr13_intensity','chr14_intensity',...
        'chr15_intensity','chr16_intensity','chr17_intensity','chr18_intensity',...
        'chr19_intensity','chrX_intensity'});
    
    % convert into nm
    d.x = d.x * xPixToNm;
    d.y = d.y * yPixToNm;
    d.z = d.z * zPixToNm;
    
    d = renamevars(d,'regionID_hyb1_60_','regionID');
end
