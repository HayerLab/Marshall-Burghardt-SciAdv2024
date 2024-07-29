%% Kinetics of RhoGTPase activities following drug treatments
%
% 170810 - Rho,Rac,Cdc42 in response to Thrombin
% 170816 - Rho,Rac,Cdc42 in response to thrombin

% Load data from all  experiments. 
% Normalize/correct for bleaching  
% plote means +/- CI

clear; close all; clc;
analysisdir='..\analysis';


%% Experiment 1. 170810 (Thrombin)
datadir=['..\170810_3i\Thrombin\data'];
sensors={
    'Rho'
    'Rac'
    'Cdc42'};
conditions={
    'Thrombin'};

for k=1:size(sensors,1)
    files=getFilenames(datadir,'_Bleach_raw');
    files=files(boolRegExp(files,sensors{k}));
    ControlFiles=files(boolRegExp(files,'control'));
    ThrombinFiles=files(boolRegExp(files,'Thrombin'));
    
    control=[];controlNorm=[];controlNormNorm=[];
    %Normalize by the mean of first 10 frames
    for i=1:size(ControlFiles,1)
        load([datadir,filesep,ControlFiles{i}]);
        control(:,i)=bleach_raw';
        controlNorm(:,i)=(bleach_raw/mean(bleach_raw(1:10)))';
      
    end
    %Normalize by mean of normalized controls
    for i=1:size(ControlFiles,1)
        controlNormNorm(:,i)=controlNorm(:,i)./mean(controlNorm,2);
        plot(controlNormNorm(:,i));hold on;
    end
    
    Thrombin=[];ThrombinNorm=[];ThrombinNormNorm=[];
    for i=1:size(ThrombinFiles,1)
        load([datadir,filesep,ThrombinFiles{i}]);
        Thrombin(:,i)=bleach_raw';
        ThrombinNorm(:,i)=(bleach_raw/mean(bleach_raw(1:10)))';
        ThrombinNormNorm(:,i)=ThrombinNorm(:,i)./mean(controlNorm,2);
        plot(ThrombinNormNorm(:,i));hold on;
    end
    
    save([analysisdir,filesep,sensors{k},'_',conditions{1},'_1.mat'],'controlNorm','controlNormNorm','ThrombinNorm','ThrombinNormNorm');
end


%% Experiment 2: 170816

datadir=['..\170816_3i\Thrombin\data'];
sensors={
    'Rho'
    'Rac'
    'Cdc42'}
conditions={
    'Thrombin'};

for k=1:size(sensors,1)
    files=getFilenames(datadir,'_Bleach_raw');
    files=files(boolRegExp(files,sensors{k}));
    ControlFiles=files(boolRegExp(files,'control'));
    ThrombinFiles=files(boolRegExp(files,'-Thrombin05'));
    
    control=[];controlNorm=[];controlNormNorm=[];
    %Normalize by the mean of first 10 frames
    for i=1:size(ControlFiles,1)
        load([datadir,filesep,ControlFiles{i}]);
        control(:,i)=bleach_raw';
        controlNorm(:,i)=(bleach_raw/mean(bleach_raw(1:10)))';
      
    end
    %Normalize by mean of normalized controls
    for i=1:size(ControlFiles,1)
        controlNormNorm(:,i)=controlNorm(:,i)./mean(controlNorm,2);
        plot(controlNormNorm(:,i));hold on;
    end
    
    Thrombin=[];ThrombinNorm=[];ThrombinNormNorm=[];
    for i=1:size(ThrombinFiles,1)
        load([datadir,filesep,ThrombinFiles{i}]);
        Thrombin(:,i)=bleach_raw';
        ThrombinNorm(:,i)=(bleach_raw/mean(bleach_raw(1:10)))';
        ThrombinNormNorm(:,i)=ThrombinNorm(:,i)./mean(controlNorm,2);
        plot(ThrombinNormNorm(:,i));hold on;
    end
    
    save([analysisdir,filesep,sensors{k},'_',conditions{1},'_2.mat'],'controlNorm','controlNormNorm','ThrombinNorm','ThrombinNormNorm');
end



%% Consolidate and plot data
% 170816 only has 30 time-points, 170810 has 60 timepoints. Plot only 30
% time points
clear;close all;
analysisdir='..\analysis';
colors={
    [0 0 0];
    [1 0 1];
    [0 1 0]};
sensors={
    'Rho'
    'Rac'
    'Cdc42'};

control_temp=[];thrombin_temp=[];
for k=1:size(sensors,1)
    load([analysisdir,filesep,sensors{k},'_Thrombin_1.mat']);
    thrombin_temp=ThrombinNormNorm(1:30,:);
    load([analysisdir,filesep,sensors{k},'_Thrombin_2.mat']);
    thrombin_temp=[thrombin_temp,ThrombinNormNorm];
    for times=1:30
        pd=fitdist(thrombin_temp(times,:)','Normal');
        thrombinStats(times,1)=pd.mu; % mean
        thrombinCI=paramci(pd);
        thrombinStats(times,2)=thrombinCI(1,1); %lower limit CI
        thrombinStats(times,3)=thrombinCI(2,1); %upper limit CI
    end
    plot(-90:10:200,thrombinStats(:,1),'color',colors{k}, 'LineWidth',1);hold on;
    plot(-90:10:200,thrombinStats(:,2),'color',colors{k},'LineWidth',0.25);
    plot(-90:10:200,thrombinStats(:,3),'color',colors{k},'LineWidth',0.25);    
    text(-50,(1.2+0.05*k),[sensors{k},' n=',num2str(size(thrombin_temp,2))],'Color',colors{k});
    
end
axis square;    
    
xlabel('time relative to Thrombin addition [sec]');
ylabel('Relative RhoGTPase FRET/CFP');
savefig('..\Fig_1f.fig');
