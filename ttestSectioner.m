function [ section1,section2,section3 ] = ttestSectioner( metric,laneIdx,congruence)
%ttestSectioner Allows to break down metric by section (1st single lane, double lane, second single lane)
%   Achieve subject-wise average of metric for each section of trial for
%   purpose of ttest (compare congruent and incongruent performance in each section)

Section1_con = [];
Section2_con = [];
Section3_con = [];
Section1_incon = [];
Section2_incon = [];
Section3_incon = [];

for i = 1:size(metric,1)
    for j = 1:size(metric,2)
        if(~isempty(metric{i,j}))
            crntCell = metric{i,j};
            if(congruence(i,j))
                Section1_con = [Section1_con;crntCell(1:laneIdx(1))]; %concat sections of trials
                Section2_con = [Section2_con;crntCell(laneIdx(1):laneIdx(2))];
                Section3_con = [Section3_con;crntCell(laneIdx(2):end)];
            else
                Section1_incon = [Section1_incon;crntCell(1:laneIdx(1))];
                Section2_incon = [Section2_incon;crntCell(laneIdx(1):laneIdx(2))];
                Section3_incon = [Section3_incon;crntCell(laneIdx(2):end)];
            end
        end
    end
    %mean sections
    section1.mean(i,1:2) = [mean(Section1_con,1) mean(Section1_incon,1)];
    section2.mean(i,1:2) = [mean(Section2_con,1) mean(Section2_incon,1)];
    section3.mean(i,1:2) = [mean(Section3_con,1) mean(Section3_incon,1)];
    %rms sections
    section1.rms(i,1:2) = [rms(Section1_con,1) mean(Section1_incon,1)];
    section2.rms(i,1:2) = [rms(Section2_con,1) mean(Section2_incon,1)];
    section3.rms(i,1:2) = [rms(Section3_con,1) mean(Section3_incon,1)];
end

end

