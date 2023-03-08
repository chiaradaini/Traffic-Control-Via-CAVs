vel=(0:5:140);
for i=1:length(vel)
    TFC_curve(i)=Simulator_function(vel(i));
end

% TFC_curve = TFC_curve*60;

figure;
plot(vel,TFC_curve);
xlim([20,140]);
ylim([26600,27800]);
ylabel('Total fuel consumption [l/h]', 'fontsize',20);
% ylabel('Average Travel Time [min]', 'fontsize',20);
xlabel('CAV velocity [km/h]','fontsize',20);

% %%%saveas('TFC trend for V between 0 and 140.png');



% pos=(2:2:50);
% for i=1:length(pos)
%     TFC_curve(i)=Simulator_function(pos(i));
% end
% 
% TFC_curve = TFC_curve*60;
% 
% figure;
% plot(pos,TFC_curve);
% % ylim([26600,27800]);
% xlim([2,50]);
% ylabel('Total fuel consumption [l/h]', 'fontsize',20);
% % ylabel('Average Travel Time [min]', 'fontsize',20);
% xlabel('CAV initial position [km]','fontsize',20);
