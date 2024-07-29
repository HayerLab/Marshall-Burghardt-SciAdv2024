%% ParallelMaskingIX.m
% This script sets up parallel processing to generate cell masks, align imgaes
% in multiple channels 
% Arnold Hayer 180622

%% Parameter setup
clear;clc;close all; 
root = 'H:\180620_IXM\180620-I';
rawdir=[root,filesep,'raw'];
bgdir=[root,filesep,'background'];
load([bgdir,filesep,'alignment parameters pX pY.mat']);
datadir=[root,filesep,'data_180622'];
if ~exist(datadir)
    mkdir(datadir);
end


k=0;
for row=5%2:3
    for col=8%2:10
        for site=1
            k=k+1;
            position{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
        end
    end
end

%% Parallel loop
for k=1:length(position)
    getMaskDataIX(position{k},bgdir,rawdir,datadir);
end
disp('done!');

