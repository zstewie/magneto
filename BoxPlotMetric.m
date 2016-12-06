function [] = BoxPlotMetric(metric, metricMean, mytitle, metricLabel, congruence, axes)
%BOXPLOTMETRIC(metric, metricMean, mytitle, metricLabel, congruence, axes)
%Make a boxplot for given metric, only to compare congruent & incongruent trials
%   leave axes empty ([]) for auto, specify [xmin xmax ymin ymax] if
%   desired. Metric mean contains overall mean for congruent and
%   incongruent trials --> metricMean = (incongruentmean,congruentmean)

NaNFiller = NaN((length(metric(congruence == 0)) - ...
    length(metric(congruence == 1))),1);
figure;
boxplot([metric(congruence == 0) ...
    [metric(congruence == 1);NaNFiller]],...
    'notch','on','colors','rb');
hold on;
scatter(2,metricMean(2),'*');
scatter(1,metricMean(1),'*');
hold off;
title(mytitle,'FontSize',12)
xlabel('Steering Mode','FontSize',12);
ylabel(metricLabel,'FontSize',12);
set(gca,'XTickLabel',{'Incongruent','Congruent'});

if(isempty(axes))
    axis auto;
else
    axis(axes);
end

end

