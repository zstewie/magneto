function [result] = profx_stats(result)

alpha = 0.05;
numMetrics = size(result.ttestMat_congruent,2);
ttestRows = {'numReversals';'numPerSecReversals';'hwrateMean';'hwaccelMean';'hwjerkMean';...
    'latVelMean';'latAccelMean';'time2hwpeak';'time2posypeak';'dist2hwpeak';'dist2posypeak';...
    'yawrateMean';'outOfBoundsTime';'latError';'hwrateRMS';'hwaccelRMS';'hwjerkRMS';...
    'latVelRMS';'latAccelRMS';'yawrateRMS';'hwratesection1mean';'hwratesection2mean';'hwratesection3mean';...
    'hwaccelsection1mean';'hwaccelsection2mean';'hwaccelsection3mean';'hwjerksection1mean';'hwjerksection2mean';...
    'hwjerksection3mean';'yawratesection1mean';'yawratesection2mean';'yawratesection3mean'};

%% Rearrange data into table
% for i = 1:size(result.participant,1) %number of participants
%     for j = 1:(2 * numMetrics)
%
%         %reshape data, separating by congruence
%         group1 = metric(find(group.var == 0));
%         group2 = metric(find(group.var == 1));
%
%     end
% end
%
% %%
%
% %get rid of NaNs
% group1(isnan(group1)) = [];
% group2(isnan(group2)) = [];
%
% % if(length(group1) > length(group2))
% %     group1 = group1(1:length(group2));
% % else
% %     group2 = group2(1:length(group1));
% % end
%
% %transpose vertical to horizontal
% group1 = group1';
% group2 = group2';
%
% %disp(['Size of group 1: ',num2str(size(group1,1)),' by ',num2str(size(group1,2))]);
% %disp(['Size of group 2: ',num2str(size(group2,1)),' by ',num2str(size(group2,2))]);
%
% %horizontal concatenation
% groups = [group1,group2];
%
% %change labels depending on passed group
% if(strcmp(group.name,'TurnDirection'))
%     labels = {'Left','Right'};
% else
%     labels = {'Incongruent','Congruent'};
% end
%
% label1 = cell(size(group1));
% label2 = cell(size(group2));
%
% %create group label cell array of same size as the corresponding data
% for i = 1:size(group1,2)
%     label1{i} = labels{1};
% end
% for i = 1:size(group2,2)
%     label2{i} = labels{2};
% end

%% Actually do stats
%last update 08.02.16

%display what groups 1 and 2 correspond to, depending on name of grouping
%var passed in
% if strcmp(group.name,'Congruence')
%     disp('Group 1: Incongruent | Group 2: Congruent');
% else
%     disp('Group 1: Left | Group 2: Right');
% end

%test for normality for group 1 and group 2 data
% figure
% normplot(group1);
% title(['Normality Plot for ' labels{1} ' trials']);
% figure
% normplot(group2);
% title(['Normality Plot for ' labels{2} ' trials']);

%test for homoscedasticity
% res1 = group1 - mean(group1);
% res2 = group2 - mean(group2);
% [result(1).h1_arch,result(1).p1_arch,result(1).fStat1_arch,result(1).crit1_arch] = ...
%     archtest(res1,'Lags',2);
% [result(1).h2_arch,result(1).p2_arch,result(1).fStat2_arch,result(1).crit2_arch] = ...
%     archtest(res2,'Lags',2);
% 
% disp('%%%%%Test for Homoscedasticity%%%%%');
% disp(result(1));

%one-way ANOVA
% labels = [label1,label2];
% [result(2).p,result(2).tbl,result(2).stats] = anova1(groups,labels,'off');
% disp('%%%%%1-way RM ANOVA%%%%%');
% disp(['p-value = ', num2str(result(2).p)]);
% disp(result(2).tbl)

%2 sample t-test for each metric
for i = 1:numMetrics
    %[result(3).h,result(3).p,result(3).ci,result(3).stats] = ttest2(group1,group2);
    [h(i,1),p(i,1),ci{i,1},stats{i,1}] = ttest(result.ttestMat(:,i),result.ttestMat(:,i+numMetrics));
    ci_low(i,1) = ci{i}(1);
    ci_hi(i,1) = ci{i}(2);
    tstat(i,1) = stats{i,1}.tstat;
    df(i,1) = stats{i,1}.df;
    sd(i,1) = stats{i,1}.sd;
end
result.stats.ttest = table(h,p,ci_low,ci_hi,tstat,df,sd,'rowNames',ttestRows);

%False discovery rate
[result.stats.ttest_pfdr,result.stats.ttest.pmask] = fdr(p,alpha);
result.stats.ttest.pfdr = ones(size(p)) * result.stats.ttest_pfdr;


end