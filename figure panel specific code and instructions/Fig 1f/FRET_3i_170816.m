%% FRET_3i_170816.m
% This script sets up parallel processing to obtain FRET data of time-lapse
% series by calling getFRETData in a parfor loop.

%% Parameter setup

clear;clc;close all; warning off;
root = '..\170816_3i\Thrombin';
rawdir=[root,filesep,'raw'];
bgdir=[root,filesep,'background'];
load([bgdir,filesep,'alignment parameters pX pY.mat']);
datadir=[root,filesep,'data'];
if ~exist(datadir)
    mkdir(datadir);
end

files=getFilenames(rawdir);
files=files(boolRegExp(files,'.tiff'));
CFPfiles=files(boolRegExp(files,'_CFP'));
%CFPfiles=CFPfiles(boolRegExp(CFPfiles,'-Rac'));
%% Parallel loop
parfor k=1:length(CFPfiles)
    getFRETData3i_kinetics(k,bgdir,rawdir,datadir);
end
disp('done!');


