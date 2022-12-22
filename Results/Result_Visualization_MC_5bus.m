clear all
close all
clc

load MC_Result
load gPC_Result

MC_obj = table2array(MC_Result(1:71,2));
gPC_obj = table2array(gPC_Result(1:71,2));
pen_level = 100 .* table2array(gPC_Result(:,1));

diff_perc = 100 * abs(gPC_obj - MC_obj) ./ MC_obj;



%%
f = figure;
f.Position = [100 100 750 500];
font = 18;
yyaxis left
stem(pen_level(1:length(gPC_obj)), gPC_obj, '+', "filled", 'LineStyle', 'none', 'color', 'k', 'LineWidth', 1)
hold on
stem(pen_level(1:length(MC_obj)), MC_obj, '.', "filled", 'LineStyle', 'none', 'color', 'k', 'LineWidth', 1)
ylim([2000 5000])
xlabel("Penetration Level [%]")
ylabel("Expected Cost [$]")


yyaxis right
% bar_graph = bar(pen_level(1:length(diff_perc)), diff_perc, 0.5, 'k');
bar_graph = bar(pen_level(1:length(diff_perc)), diff_perc, "LineWidth",1);
bar_graph.FaceColor = 'none';
hatchfill2(bar_graph,'single','HatchAngle',45)
ylabel("Difference in Expected Costs [%]")
ylim([0 6])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gca, 'FontSize', font)



legend('Proposed Method', 'Monte-Carlo Sim.', 'Percentage Diff.')





