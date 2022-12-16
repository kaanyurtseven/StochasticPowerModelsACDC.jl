clear all
close all
clc

load Results4

Case1_obj = table2array(Results4(:,3));
Case2_obj = table2array(Results4(:,4));
pen_level = 100 .* table2array(Results4(:,1));

Case1_obj(49:end) = [];
Case1_obj(49:54) = NaN;

Case2_obj(55:end) = [];
Case2_obj(40) = 230009.2618;
Case2_obj(41) = 229239.4826;




diff = abs(Case2_obj - Case1_obj);

%%
figure;
font = 18;

plot(pen_level(1:length(Case1_obj)), Case1_obj, "LineWidth", 2)
hold on
plot(pen_level(1:length(Case2_obj)), Case2_obj, "LineWidth", 2)
xlabel("PV Hosting Level [%]")
ylabel("Expected Cost [$]")
legend('AC', 'ACDC')
set(gca, 'FontSize', font)

%location of subpart on figure
xstart=.17;
xend=.45;
ystart=.25;
yend=.45;
axes('position',[xstart ystart xend-xstart yend-ystart ])
box on
range = 45:54;
plot(pen_level(range), Case1_obj(range), 'LineWidth', 2)
hold on
plot(pen_level(range), Case2_obj(range), 'LineWidth', 2)
set(gca,'Xtick',45:1:55, 'Ytick', [], 'FontSize', font-6)
%%
figure;
plot(pen_level(1:length(diff(1:end-1))), diff(1:end-1), "LineWidth", 2)
xlabel("PV Hosting Level [%]")
ylabel("Difference in Expected Costs [$]")
set(gca, 'FontSize', font)