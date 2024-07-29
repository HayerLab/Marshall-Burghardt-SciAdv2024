%% Kinetics of RhoGTPase activities following rapamycin-incuded translocation of 
% the RhoGEFs (ARHGEF1, TIAM1, ITSN1) and the RhoGAP ARHGAP29. Combine
% replicates from experiments: 
% 170801 - Rho,Rac,Cdc42 activation
% 170810 - ARHGAP29 activation
% 170920 - Rho, Rac, Cdc42, ARHGAP29 activation

% N.B: There was a labeling error in 170920: RacAct is actually RhoAct and RhoAct is RacAct.

% Load data from all  experiments. Normalize by control cells/correct for bleaching  

clear; close all; clc;
analysisdir='..\analysis';

%% Experiment 1: 170801
datadir='..\170801_3i\data';
sensors={
    'Rho_' 
    'Rac_'
    'Cdc42_'};
conditions={
    '1087_'     % Rho activation
    '1103_'     % Rac activation
    '1154_'};   % Cdc42 activation
labels={
    'RhoAct'
    'RacAct'
    'Cdc42Act'};
for k=1:size(sensors,1)
    for n=1:size(conditions,1)
        files=getFilenames(datadir,'_kinetics');
        files=files(boolRegExp(files,sensors{k}));
        files=files(boolRegExp(files,conditions{n}));
        ROIcontrolAll=[];ROItestAll=[];ROIcontrolAllNorm=[];ROItestAllNorm=[];
        for j=1:size(files,1)
            load([datadir,filesep,files{j}]);
            ROIcontrolAll=[ROIcontrolAll,ROIcontrol];
            ROItestAll=[ROItestAll,ROItest];
        end
        % normalize each trace by the mean of first 10 timepoints
        for i=1:size(ROIcontrolAll,2)
            ROIcontrolAllNorm(:,i)=ROIcontrolAll(:,i)/mean(ROIcontrolAll(1:10,i));
        end
        ROIcontrolMean=mean(ROIcontrolAllNorm,2);
        for i=1:size(ROIcontrolAll,2)
            ROIcontrolAllNorm(:,i)=ROIcontrolAllNorm(:,i)./ROIcontrolMean;
        end
        for i=1:size(ROItestAll,2)
            ROItestAllNorm(:,i)=ROItestAll(:,i)/mean(ROItestAll(1:10,i));
            ROItestAllNorm(:,i)=ROItestAllNorm(:,i)./ROIcontrolMean;
        end
        save([analysisdir,filesep,sensors{k}(1:end-1),labels{n},'_1.mat'],...
               'ROIcontrolAll','ROItestAll','ROIcontrolAllNorm','ROItestAllNorm');
    end
end
%% Experiment 170810 (ARHAGP29 only)
datadir='..\170810_3i\data';
sensors={
    'Rho'
    'Rac'
    'Cdc42'};
conditions={
    'GAP29WT'};
labels={
    'AG29Act'};
for k=1:size(sensors,1)
    for n=1:size(conditions,1)
        files=getFilenames(datadir,'_kinetics');
        files=files(boolRegExp(files,sensors{k}));
        files=files(boolRegExp(files,conditions{n}));
        ROIcontrolAll=[];ROItestAll=[];ROIcontrolAllNorm=[];ROItestAllNorm=[];
        for j=1:size(files,1)
            load([datadir,filesep,files{j}]);
            ROIcontrolAll=[ROIcontrolAll,ROIcontrol];
            ROItestAll=[ROItestAll,ROItest];
        end
        % normalize each trace by the mean of first 10 timepoints
        for i=1:size(ROIcontrolAll,2)
            ROIcontrolAllNorm(:,i)=ROIcontrolAll(:,i)/mean(ROIcontrolAll(1:10,i));
        end
        ROIcontrolMean=mean(ROIcontrolAllNorm,2);
        for i=1:size(ROIcontrolAll,2)
            ROIcontrolAllNorm(:,i)=ROIcontrolAllNorm(:,i)./ROIcontrolMean;
        end
        for i=1:size(ROItestAll,2)
            ROItestAllNorm(:,i)=ROItestAll(:,i)/mean(ROItestAll(1:10,i));
            ROItestAllNorm(:,i)=ROItestAllNorm(:,i)./ROIcontrolMean;
        end
        save([analysisdir,filesep,sensors{k},labels{n},'_1.mat'],...
               'ROIcontrolAll','ROItestAll','ROIcontrolAllNorm','ROItestAllNorm');
    end
end
%% Experiment 2: 170920
datadir='..\170920_3i\data';
sensors={
    'Rho_'
    'Rac_'
    'Cdc42_'};
conditions={
    'RhoAct'
    'RacAct'
    'Cdc42Act'
    'AG29Act'};
labels={
    'RacAct' % Rho/rac swapped due to labeling error
    'RhoAct' % Rho/rac swapped due to labeling error
    'Cdc42Act'
    'AG29Act'};
for k=1:size(sensors,1)
    for n=1:size(conditions,1)
        files=getFilenames(datadir,'_kinetics');
        files=files(boolRegExp(files,sensors{k}));
        files=files(boolRegExp(files,conditions{n}));
        ROIcontrolAll=[];ROItestAll=[];ROIcontrolAllNorm=[];ROItestAllNorm=[];
        for j=1:size(files,1)
            load([datadir,filesep,files{j}]);
            ROIcontrolAll=[ROIcontrolAll,ROIcontrol];
            ROItestAll=[ROItestAll,ROItest];
        end
        % normalize each trace by the mean of first 10 timepoints
        for i=1:size(ROIcontrolAll,2)
            ROIcontrolAllNorm(:,i)=ROIcontrolAll(:,i)/mean(ROIcontrolAll(1:10,i));
        end
        ROIcontrolMean=mean(ROIcontrolAllNorm,2);
        for i=1:size(ROIcontrolAll,2)
            ROIcontrolAllNorm(:,i)=ROIcontrolAllNorm(:,i)./ROIcontrolMean;
        end
        for i=1:size(ROItestAll,2)
            ROItestAllNorm(:,i)=ROItestAll(:,i)/mean(ROItestAll(1:10,i));
            ROItestAllNorm(:,i)=ROItestAllNorm(:,i)./ROIcontrolMean;
        end
        save([analysisdir,filesep,sensors{k}(1:end-1),labels{n},'_2.mat'],...
               'ROIcontrolAll','ROItestAll','ROIcontrolAllNorm','ROItestAllNorm');
    end
end
%% Consolidate data
clear;
analysisdir='..\analysis';
conditions={
    'RhoRhoAct'
    'RacRhoAct'
    'Cdc42RhoAct'
    'RhoRacAct'
    'RacRacAct'
    'Cdc42RacAct'
    'RhoCdc42Act'
    'RacCdc42Act'
    'Cdc42Cdc42Act'
    'RhoAG29Act'
    'RacAG29Act'
    'Cdc42AG29Act'
    };
 
for k=1:size(conditions,1)
    control_temp=[];test_temp=[];
    load([analysisdir,filesep,conditions{k},'_1.mat']);
    control_temp=ROIcontrolAllNorm;
    test_temp=ROItestAllNorm;
    load([analysisdir,filesep,conditions{k},'_2.mat']);
    control_temp=[control_temp,ROIcontrolAllNorm];
    test_temp=[test_temp,ROItestAllNorm];

    subplot(4,3,k);
  
    plot(-135:15:600,mean(control_temp,2),'color',[0 0 0], 'LineWidt',2);hold on;
    % With +/-SD
    %plot(-135:15:600,mean(control_temp,2)+std(control_temp,0,2),'color',[0 0 0], 'LineWidt',1);hold on;
    %plot(-135:15:600,mean(control_temp,2)-std(control_temp,0,2),'color',[0 0 0], 'LineWidt',1);hold on;
    % With +/-SEM:
    %plot(-135:15:600,mean(control_temp,2)+std(control_temp,0,2)/sqrt(size(control_temp,2)),'color',[0 0 0], 'LineWidt',1);hold on;
    %plot(-135:15:600,mean(control_temp,2)-std(control_temp,0,2)/sqrt(size(control_temp,2)),'color',[0 0 0], 'LineWidt',1);hold on;
    % With 95% CI:
    plot(-135:15:600,mean(control_temp,2)+1.96*std(control_temp,0,2)/sqrt(length(control_temp)),'color',[0 0 0], 'LineWidt',1);hold on;
    plot(-135:15:600,mean(control_temp,2)-1.96*std(control_temp,0,2)/sqrt(length(control_temp)),'color',[0 0 0], 'LineWidt',1);hold on;
    
    
    plot(-135:15:600,mean(test_temp,2),'color',[1 0 0], 'LineWidt',2);
        
    % With +/- SD:
    %  plot(-135:15:600,mean(test_temp,2)+std(test_temp,0,2),'color',[1 0 0], 'LineWidt',1);hold on;
    %  plot(-135:15:600,mean(test_temp,2)-std(test_temp,0,2),'color',[1 0 0], 'LineWidt',1);hold on;
    % With +/- SEM:
    %  plot(-135:15:600,mean(test_temp,2)+std(test_temp,0,2)/sqrt(size(test_temp,2)),'color',[1 0 0], 'LineWidt',1);hold on;
    %  plot(-135:15:600,mean(test_temp,2)-std(test_temp,0,2)/sqrt(size(test_temp,2)),'color',[1 0 0], 'LineWidt',1);hold on;
    % With 95% CI:
    plot(-135:15:600,mean(test_temp,2)+1.96*std(test_temp,0,2)/sqrt(length(test_temp)),'color',[1 0 0], 'LineWidt',1);hold on;
    plot(-135:15:600,mean(test_temp,2)-1.96*std(test_temp,0,2)/sqrt(length(test_temp)),'color',[1 0 0], 'LineWidt',1);hold on;
    
    
    hold off;
    xlim([-150 600]);ylim([0.95 1.25]);title(conditions{k});
    xlabel(['n=',num2str(size(control_temp,2)),'(control)  ','n=',num2str(size(test_temp,2)),'(test)']);
    
    title(conditions{k});
    
end


%% Plot data;
close all;
% Rho activation, plot Rho, Rac, Cdc42 kinetics. 
xrange=[-150 600];
yrange=[0.8 1.4];
colors={
    [0 0 0];
    [1 0 1];
    [0 1 0]};

subplot(2,2,1);
conditions={
    'RhoRhoAct'
    'RacRhoAct'
    'Cdc42RhoAct'};
test_tempStats=[];test_temp=[];
for k=1:length(conditions)
    load([analysisdir,filesep,conditions{k},'_1.mat']);
    test_temp=ROItestAllNorm;
    load([analysisdir,filesep,conditions{k},'_2.mat']);
    test_temp=[test_temp,ROItestAllNorm];
        
%     % With +/- SD:
%     plot(-135:15:600,mean(test_temp,2),'color',colors{k}, 'LineWidth',2);hold on;
%     plot(-135:15:600,mean(test_temp,2)+std(test_temp,0,2),'--','color',colors{k},'LineWidth',0.25);
%     plot(-135:15:600,mean(test_temp,2)-std(test_temp,0,2),'--','color',colors{k},'LineWidth',0.25);
    % With +/- CI
    % Compute 95% CI
    for times=1:50
        pd=fitdist(test_temp(times,:)','Normal');
        test_tempStats(times,1)=pd.mu; % mean
        test_tempCI=paramci(pd);
        test_tempStats(times,2)=test_tempCI(1,1); %lower limit CI
        test_tempStats(times,3)=test_tempCI(2,1); %upper limit CI
    end
    plot(-135:15:600,test_tempStats(:,1),'color',colors{k}, 'LineWidth',1);hold on;
    plot(-135:15:600,test_tempStats(:,2),'color',colors{k},'LineWidth',0.25);
    plot(-135:15:600,test_tempStats(:,3),'color',colors{k},'LineWidth',0.25);    
    text(0,(0.8+0.05*k),[conditions{k},' n=',num2str(size(test_temp,2))],'Color',colors{k});
    axis square;
end
hold off;
%xlim(xrange);ylim(yrange);title('Rho activation');
xlim(xrange);ylim(yrange);title('Rho activation');



% Rho inactivation activation, plot Rho, Rac, Cdc42 kinetics in one graph. 
subplot(2,2,2);
conditions={
    'RhoAG29Act'
    'RacAG29Act'
    'Cdc42AG29Act'};
test_tempStats=[];test_temp;
for k=1:length(conditions)
    load([analysisdir,filesep,conditions{k},'_1.mat']);
    test_temp=ROItestAllNorm;
    load([analysisdir,filesep,conditions{k},'_2.mat']);
    test_temp=[test_temp,ROItestAllNorm];
    
    % With +/- SD:
%     plot(-135:15:600,mean(test_temp,2),'color',colors{k}, 'LineWidth',2);hold on;
%     plot(-135:15:600,mean(test_temp,2)+std(test_temp,0,2),'--','color',colors{k}, 'LineWidth',0.25);
%     plot(-135:15:600,mean(test_temp,2)-std(test_temp,0,2),'--','color',colors{k}, 'LineWidth',0.25);
     % Compute 95% CI
    for times=1:50
        pd=fitdist(test_temp(times,:)','Normal');
        test_tempStats(times,1)=pd.mu; % mean
        test_tempCI=paramci(pd);
        test_tempStats(times,2)=test_tempCI(1,1); %lower limit CI
        test_tempStats(times,3)=test_tempCI(2,1); %upper limit CI
    end
    plot(-135:15:600,test_tempStats(:,1),'color',colors{k}, 'LineWidth',1);hold on;
    plot(-135:15:600,test_tempStats(:,2),'color',colors{k},'LineWidth',0.25);
    plot(-135:15:600,test_tempStats(:,3),'color',colors{k},'LineWidth',0.25);    
    text(0,(1.2+0.05*k),[conditions{k},' n=',num2str(size(test_temp,2))],'Color',colors{k});
    axis square;

end
xlim(xrange);ylim(yrange);title('ARHGAP29 activation');

% Rac activation, plot Rho, Rac, Cdc42 kinetics in one graph. 
subplot(2,2,3);
conditions={
    'RhoRacAct'
    'RacRacAct'
    'Cdc42RacAct'};
test_tempStats=[];test_temp;
for k=1:length(conditions)
    load([analysisdir,filesep,conditions{k},'_1.mat']);
    test_temp=ROItestAllNorm;
    load([analysisdir,filesep,conditions{k},'_2.mat']);
    test_temp=[test_temp,ROItestAllNorm];
    
    % With +/- SD:
%     plot(-135:15:600,mean(test_temp,2),'color',colors{k}, 'LineWidth',2);hold on;
%     plot(-135:15:600,mean(test_temp,2)+std(test_temp,0,2),'--','color',colors{k}, 'LineWidth',0.25);
%     plot(-135:15:600,mean(test_temp,2)-std(test_temp,0,2),'--','color',colors{k}, 'LineWidth',0.25);
     % Compute 95% CI
    for times=1:50
        pd=fitdist(test_temp(times,:)','Normal');
        test_tempStats(times,1)=pd.mu; % mean
        test_tempCI=paramci(pd);
        test_tempStats(times,2)=test_tempCI(1,1); %lower limit CI
        test_tempStats(times,3)=test_tempCI(2,1); %upper limit CI
    end
    plot(-135:15:600,test_tempStats(:,1),'color',colors{k}, 'LineWidth',1);hold on;
    plot(-135:15:600,test_tempStats(:,2),'color',colors{k},'LineWidth',0.25);
    plot(-135:15:600,test_tempStats(:,3),'color',colors{k},'LineWidth',0.25);    
    text(0,(0.8+0.05*k),[conditions{k},' n=',num2str(size(test_temp,2))],'Color',colors{k});
    axis square;
end
xlim(xrange);ylim(yrange);title('Rac activation');
 
% Cdc42 activation, plot Rho, Rac, Cdc42 kinetics in one graph. 
subplot(2,2,4);
conditions={
    'RhoCdc42Act'
    'RacRacAct'
    'Cdc42Cdc42Act'};
test_tempStats=[];test_temp;
for k=1:length(conditions)
    load([analysisdir,filesep,conditions{k},'_1.mat']);
    test_temp=ROItestAllNorm;
    load([analysisdir,filesep,conditions{k},'_2.mat']);
    test_temp=[test_temp,ROItestAllNorm];
    
    % With +/- SD:
%     plot(-135:15:600,mean(test_temp,2),'color',colors{k}, 'LineWidth',2);hold on;
%     plot(-135:15:600,mean(test_temp,2)+std(test_temp,0,2),'--','color',colors{k}, 'LineWidth',0.25);
%     plot(-135:15:600,mean(test_temp,2)-std(test_temp,0,2),'--','color',colors{k}, 'LineWidth',0.25);
     % Compute 95% CI
    for times=1:50
        pd=fitdist(test_temp(times,:)','Normal');
        test_tempStats(times,1)=pd.mu; % mean
        test_tempCI=paramci(pd);
        test_tempStats(times,2)=test_tempCI(1,1); %lower limit CI
        test_tempStats(times,3)=test_tempCI(2,1); %upper limit CI
    end
    plot(-135:15:600,test_tempStats(:,1),'color',colors{k}, 'LineWidth',1);hold on;
    plot(-135:15:600,test_tempStats(:,2),'color',colors{k},'LineWidth',0.25);
    plot(-135:15:600,test_tempStats(:,3),'color',colors{k},'LineWidth',0.25);    
    text(0,(0.8+0.05*k),[conditions{k},' n=',num2str(size(test_temp,2))],'Color',colors{k});
    axis square;
end
xlim(xrange);ylim(yrange);title('Cdc42 activation');












