function [ttestGroup1,ttestGroup2] = ttestFormat(metric,separator,separatorName)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dim1 = size(metric,1);
dim2 = size(metric,2);
separator = separator(1:dim1,1:dim2);

if(strcmp(separatorName,'congruence'))
    for i = 1:size(metric,1)
        crntRow = metric(i,:);
        crntIncongruent = crntRow(separator(i,:) == 0);
        crntCongruent = crntRow(separator(i,:) == 1);
        IncongruentMean(i,1) = mean(crntIncongruent,'omitnan');
        CongruentMean(i,1) = mean(crntCongruent,'omitnan');
    end
    ttestGroup1 = IncongruentMean;
    ttestGroup2 = CongruentMean;
elseif(strcmp(separatorName,'group'))
    for i = 1:size(metric,1)
        crntRow = metric(i,:);
        crntGroup1 = crntRow(separator(i,:) == 1);
        crntGroup2 = crntRow(separator(i,:) == 2);
        Group1Mean(i,1) = mean(crntGroup1,'omitnan');
        Group2Mean(i,1) = mean(crntGroup2,'omitnan');
    end
    ttestGroup1 = Group1Mean;
    ttestGroup2 = Group2Mean;
end
end

