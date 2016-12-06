function [ metricOut ] = metricSplitter( metricIn,trialType )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

group1 = (trialType == 1);
group2 = (trialType == 2);

metricOut.group1 = NaN(size(trialType));
metricOut.group2 = NaN(size(trialType));

%signal
metricOut.group1 = metricIn(group1);
metricOut.group2 = metricIn(group2);

%mean
metricOut.group1Mean = mean(metricOut.group1,'omitnan');
metricOut.group2Mean = mean(metricOut.group2,'omitnan');

%standard deviation
metricOut.group1Std = std(metricOut.group1,'omitnan');
metricOut.group2Std = std(metricOut.group2,'omitnan');

%rms
metricOut.group1RMS = rms(metricOut.group1,'omitnan');
metricOut.group2RMS = rms(metricOut.group2,'omitnan');

end

