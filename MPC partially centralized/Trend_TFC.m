%Results with 6 minuts of optimization and 5 of simulation CAV
TFC_duplex_centralized_MPC_6min = 1.0e+04 * [2.713192757346654   2.700751693859717   2.689634590818775   2.678784275198190   2.670653864810229 2.662969004411471   2.658574856556342   2.655115330898539   2.651489932404500   2.652732093959461];
TFC_duplex_decentralized_MPC_6min = 1.0e+04 * [2.712629796261012   2.701020281871168   2.697586730951433   2.683585284712902   2.685824791489196 2.677382725999736   2.674949672049186   2.671868109307289   2.665360183433344   2.663873170446402];
TFC_duplex_quasidecentralized_MPC_6min = 1.0e+04 * [2.712651194293895   2.702283199093428   2.697962510579064   2.681933242447742   2.678115501993395 2.676086541130631   2.666161297434678   2.666975778927724   2.665026717763898   2.666161424873439];

%Results with 15 minuts of optimization and 5 of simulation CAV
TFC_duplex_centralized_MPC_15min = 1.0e+04 * [2.695168004285074   2.686657777345800   2.680105349842892   2.673372154551031   2.670603773913582 2.668312057550049   2.657038930312783   2.654287840643390   2.651951259298137   2.645900279869015];
TFC_duplex_decentralized_MPC_15min = 1.0e+04 * [2.695187956775711   2.696601521016823   2.690782889816657   2.686152430288735   2.682983729302322 2.676511304519656   2.674287802711610   2.670225003574413   2.666536300932551   2.668223991938429];
TFC_duplex_quasidecentralized_MPC_15min = 1.0e+04 * [2.696341379762262   2.689762228996084   2.688910441377225   2.684468122092172   2.679812062189566   2.668273325842654   2.661481840980211   2.658289160054456   2.655324443590476   2.652714310924141];

x = 1:10;

TFC_MPC_poly_cent_6min = fit(x.',TFC_duplex_centralized_MPC_6min.','poly2');
TFC_MPC_poly_dec_6min = fit(x.',TFC_duplex_decentralized_MPC_6min.','poly2');
TFC_MPC_poly_quasidec_6min = fit(x.',TFC_duplex_quasidecentralized_MPC_6min.','poly2');

TFC_MPC_poly_cent_15 = fit(x.',TFC_duplex_centralized_MPC_15min.','poly2');
TFC_MPC_poly_dec_15 = fit(x.',TFC_duplex_decentralized_MPC_15min.','poly2');
TFC_MPC_poly_quasidec_15 = fit(x.',TFC_duplex_quasidecentralized_MPC_15min.','poly2');

hold on
figure(1);
TFCnoMPC = (2.7647*1e4)*ones(1,length(TFC_duplex_centralized_MPC_6min));
plot(TFCnoMPC, 'k', 'linewidth', 2);

plot(TFC_duplex_centralized_MPC_6min,'or');
plot(TFC_duplex_decentralized_MPC_6min,'ob');
plot(TFC_duplex_quasidecentralized_MPC_6min,'og');

plot(TFC_MPC_poly_cent_6min,'r');
plot(TFC_MPC_poly_dec_6min,'b');
plot(TFC_MPC_poly_quasidec_6min,'g');

% plot(TFC_duplex_centralized_MPC_15min,'or');
% plot(TFC_duplex_decentralized_MPC_15min,'ob');
% plot(TFC_duplex_quasidecentralized_MPC_15min,'og');
% 
% plot(TFC_MPC_poly_cent_15,'r');
% plot(TFC_MPC_poly_dec_15,'b');
% plot(TFC_MPC_poly_quasidec_15,'g');

% title('TFC optimisations varying the number of vehicles','fontsize',14);
xlabel('Number of controlled CAVs','fontsize',15);
ylabel('Total fuel consumption [l/h]', 'fontsize',15);
xlim([1 10]);
% legend('TFC without MPC', 'TFC 1 hour optim. centralized', 'TFC with MPC optim. centralized','TFC 1 hour optim. decentralized', 'TFC with MPC optim. decentralized', 'TFC 1 hour optim. partially centralized', 'TFC with MPC optim. partially centralized');
%legend('TFC without MPC', 'TFC MPC centralized','TFC MPC decentralized','TFC MPC partially centralized', 'fontsize',15);
hold off

% % % % saveas(gca,'TFC_optimisations.png');
