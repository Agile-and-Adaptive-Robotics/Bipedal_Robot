% clear; clc; close all;
load data20mm.mat

% Assuming x, y, z are column vectors of the same length
n = length(Ax);
rng('default'); % For reproducibility
idx = randperm(n);

% 80% training, 20% validation
numTrain = round(0.8 * n);
trainIdx = idx(1:numTrain);
valIdx = idx(numTrain+1:end);

% Training data
xTrain20 = Ax(trainIdx);
yTrain20 = Ay(trainIdx);
zTrain20 = Az(trainIdx);

% Validation data
xVal20 = Ax(valIdx);
yVal20 = Ay(valIdx);
zVal20 = Az(valIdx);

save data20mm_sorted.mat Ax Ay Az...
                        xTrain20 yTrain20 zTrain20...
                        xVal20 yVal20 zVal20