clc
close all
clear all

% physicochemical and aac matrices
load('..\Data\Variabili\Data.mat')
% load the precalculated DCT sequences for the train set
load('..\Data\Variabili\1\trainDCT.mat')
% load the precalculated DCT sequences for the test set
load('..\Data\Variabili\1\testDCT.mat')
% save the dayhoff matrix for the pssm function
dayhoffmatrix = dayhoff;


%load AMP dataset
datasetA = fastaread("..\Data\dataset\pos710.fasta");
% load nonAMP dataset
datasetNA = fastaread("..\Data\dataset\neg710.fasta");
% create a vector of train labels
trainlabel = zeros(length(datasetA)*2,1);

% variable for get the same element of the dataset
index = 1;
% fill a cell array with the traduced element of the datasets
for i = 1:2:length(datasetA)+length(datasetNA)
   % insert AMP
   Xtrain(i,1) = {aa2int(datasetA(index).Sequence)};
   % AMP label = 1
   trainlabel(i) = 1;
   % insert NONAMP
   Xtrain(i+1,1) = {aa2int(datasetNA(index).Sequence)};
   % NONAMP label = 2
   trainlabel(i+1) = 2;
   % next element
   index = index + 1;
end
clear index datasetNA datasetA

% compute the frequency matrix for each sequence of Xtrain
for i = 1:length(Xtrain)
   freq(i,1) = {Extraction.freq_vector(Xtrain{i,1}(1,:))};
end

for i = 1:length(Xtrain)
   % compute AAC, PSSM and PP matrix
   AAC(i,1) = {Extraction.aac(freq{i,1}(:,:),contact)};
   PSSM(i,1) = {Extraction.pssm(freq{i,1}(:,:),dayhoffmatrix)};
   PP(i,1) = {Extraction.physicochem(freq{i,1}(:,:),physicochemical)};
   % get the PseAAC vector using grey model
   PSeAA(i,1) = {Extraction.pseaac(Xtrain{i,1},physicochemical)};
end

% apply each extraction method at the AAC, PSSM and PP matrices
for i = 1:length(Xtrain)
   %precalculated data because the computation require more time
   %DCT_AAC(i,1) = {Extraction.discrete_wavelet(AAC{i,1}(:,:))}; % giusto
   %DCT_PSSM(i,1) = {Extraction.discrete_wavelet(PSSM{i,1}(:,:))}; % giusto
   %DCT__PP(i,1) = {Extraction.discrete_wavelet(PP{i,1}(:,:))}; % giusto
   PSE_AAC(i,1) = {Extraction.pseudo_alg(AAC{i,1}(:,:))}; 
   PSE_PSSM(i,1) = {Extraction.pseudo_alg(PSSM{i,1}(:,:))};
   PSE__PP(i,1) = {Extraction.pseudo_alg(PP{i,1}(:,:))};
   AVB_AAC(i,1) = {Extraction.AvBlock(AAC{i,1}(:,:))};
   AVB_PSSM(i,1) = {Extraction.AvBlock(PSSM{i,1}(:,:))};
   AVB__PP(i,1) = {Extraction.AvBlock(PP{i,1}(:,:))};
end

% insert all the sequence in a table
XTrain = table(trainDCT_AAC,trainDCT_PSSM,trainDCT__PP,PSE_AAC,PSE_PSSM,PSE__PP,AVB_AAC,AVB_PSSM,AVB__PP,PSeAA,trainlabel, ....
    'VariableNames', {'DCT_AAC','DCT_PSSM','DCT__PP','PSE_AAC','PSE_PSSM','PSE__PP','AVB_AAC','AVB_PSSM','AVB__PP', 'PSeAA', 'ClassNames'});
clear trainDCT_AAC trainDCT_PSSM trainDCT__PP PSE_AAC PSE_PSSM PSE__PP AVB_AAC AVB_PSSM AVB__PP PSeAA AAC PSSM PP freq Xtrain

% repeat the same thing for the test dataset

%load AMP dataset
datasetA = fastaread("Data\dataset\amp920.fasta");
% load nonAMP dataset
datasetNA = fastaread("Data\dataset\nonamp920.fasta");
% create a vector of test labels
testlabel = zeros(length(datasetA)*2,1);

index = 1;
% fill a cell array with the traduced element of the datasets
for i = 1:2:length(datasetA)+length(datasetNA)
   % insert AMP
   Xtrain(i,1) = {aa2int(datasetA(index).Sequence)};
   % AMP label = 2
   testlabel(i) = 1;
   % insert NONAMP
   Xtrain(i+1,1) = {aa2int(datasetNA(index).Sequence)};
   % NONAMP label = 1
   testlabel(i+1) = 2;
   % next element
   index = index + 1;
end
clear index datasetNA datasetA

% compute the frequency matrix for each sequence of Xtest
for i = 1:length(Xtrain)
   freq(i,1) = {Extraction.freq_vector(Xtrain{i,1}(1,:))};
end

for i = 1:length(Xtrain)
   % compute AAC, PSSM and PP matrix
   AAC(i,1) = {Extraction.aac(freq{i,1}(:,:),contact)};
   PSSM(i,1) = {Extraction.pssm(freq{i,1}(:,:),dayhoffmatrix)};
   PP(i,1) = {Extraction.physicochem(freq{i,1}(:,:),physicochemical)};
   PSeAA(i,1) = {Extraction.pseaac(Xtrain{i,1},physicochemical)};
end

% apply each extraction method at the AAC, PSSM and PP matrices
for i = 1:length(Xtrain)
   %DCT_AAC(i,1) = {Extraction.discrete_wavelet(AAC{i,1}(:,:))}; % giusto
   %DCT_PSSM(i,1) = {Extraction.discrete_wavelet(PSSM{i,1}(:,:))}; % giusto
   %DCT__PP(i,1) = {Extraction.discrete_wavelet(PP{i,1}(:,:))}; % giusto
   PSE_AAC(i,1) = {Extraction.pseudo_alg(AAC{i,1}(:,:))}; 
   PSE_PSSM(i,1) = {Extraction.pseudo_alg(PSSM{i,1}(:,:))};
   PSE__PP(i,1) = {Extraction.pseudo_alg(PP{i,1}(:,:))};
   AVB_AAC(i,1) = {Extraction.AvBlock(AAC{i,1}(:,:))};
   AVB_PSSM(i,1) = {Extraction.AvBlock(PSSM{i,1}(:,:))};
   AVB__PP(i,1) = {Extraction.AvBlock(PP{i,1}(:,:))};
end

% insert all the sequence in a table
XTest = table(testDCT_AAC,testDCT_PSSM,testDCT__PP,PSE_AAC,PSE_PSSM,PSE__PP,AVB_AAC,AVB_PSSM,AVB__PP,PSeAA,testlabel, ....
    'VariableNames', {'DCT_AAC','DCT_PSSM','DCT__PP','PSE_AAC','PSE_PSSM','PSE__PP','AVB_AAC','AVB_PSSM','AVB__PP', 'PSeAA', 'ClassNames'});
clear trainDCT_AAC trainDCT_PSSM trainDCT__PP PSE_AAC PSE_PSSM PSE__PP AVB_AAC AVB_PSSM AVB__PP PSeAA AAC PSSM PP freq Xtrain

% set net options
metodoOptim='adam';
miniBatchSize = 30;
options = trainingOptions(metodoOptim,...
    'MiniBatchSize',miniBatchSize,...
    'MaxEpochs',10,...
    'Verbose',false,...
    'Plots','training-progress');
clear metodoOptim miniBatchSize

% create layer without the sequenceInputLayer. It will be insert later
lay = [ ...
    flattenLayer("Name","flatten")
    bilstmLayer(100,'OutputMode','last',"Name","lstm")
    fullyConnectedLayer(2,"Name","conn")
    softmaxLayer("Name","softmax")
    classificationLayer("Name","class")];

% save all the score of each network
saveScore = [];
for i = 1:size(XTrain,2)-1
    % insert sequence input layer with the right input dimension for each
    % set of data
    layers = [ ...
        sequenceInputLayer([size(XTrain{1,i}{1,1}) 1], "Name", "seq")
        lay
    ];
   % extract element for training
   Xtr = table2cell(XTrain(:,i));
   % extract element for testing
   Xte = table2cell(XTest(:,i));
   % train the network
   net(i) = trainNetwork(Xtr, categorical(trainlabel), layerGraph(layers),options);
   
   [Ypred, score] = classify(net(i), Xte);
   % save all the prediction of each network
   savePred(:,i) = Ypred(:,1);
   
   % save the scores of each net. The prediction of the class 1 of net(i)
   % is in the i-th column, the prediction of the class 2 of net (i) is in
   % the (i+10)-th column
   saveScore(:,i) = score(:,1);
   saveScore(:,size(XTrain,2) + i) = score(:,2);
   
   % create confusion matrix
    cm = confusionmat(Ypred,categorical(testlabel));
    % extract True Positive, True Negative, False Positive, False Negative
    TP = cm(1,1);
    TN = cm(2,2);
    FP = cm(2,1);
    FN = cm(1,2);
    
    % calculation of sensibility, specificity, accuracy and Matthews 
    % correlation coefficient
    SENS(i) = (TP / (TP + FN));
    SPEC(i) = (TN / (TN + FP));
    ACC(i) = (TP + TN) / (TP + FP + TN + FN);
    MCC(i) = ((TP * TN) - (FN * FP)) / sqrt((TP + FN) * (TN + FP) * (TP + FP) * (TN + FN));
end

% combined prediction of all networks
finalPred = zeros([length(testlabel) 1]);
% ammply major vote rule
for i = 1:length(testlabel)
    % counter for each class
   class1 = 0;
   class2 = 0;
   for k = 1:length(ACC)
      % check which score of the k-th method is greater and unpdate the
      % counter
      if saveScore(i,k) > saveScore(i,k+size(XTrain,2))
         class1 = class1 + 1; 
      else
         class2 = class2 + 1; 
      end
   end
   % assign the element at the winning class
   if class1 > class2
       finalPred(i) = 1;
   else
       finalPred(i) = 2;
   end
end

clear ACC MCC SENS SPEC

% cast the predict to categorical
finalPred = categorical(finalPred);

cm = confusionmat(finalPred,categorical(testlabel));
% extract True Positive, True Negative, False Positive, False Negative
TP = cm(1,1);
TN = cm(2,2);
FP = cm(2,1);
FN = cm(1,2);

% calculation of sensibility, specificity, accuracy and Matthews 
% correlation coefficient
fSENS = (TP / (TP + FN));
fSPEC = (TN / (TN + FP));
fACC = (TP + TN) / (TP + FP + TN + FN);
fMCC = ((TP * TN) - (FN * FP)) / sqrt((TP + FN) * (TN + FP) * (TP + FP) * (TN + FN));

plotconfusion(categorical(testlabel),finalPred);

Qstatic = Extraction.qstatistic(savePred,testlabel);

disp(fACC);