%Results with 6 minuts of optimization and 5 of simulation
TFC_duplex_partiallycentralized_MPC_5 = 1.0e+04 * [2.712790527042111 2.701516280502553 2.692261072228699 2.686426455520727 2.677881272630511 2.673501978883587 2.673479739193135 2.663341723347655 2.664041979988226 2.667003735902267];
TFC_duplex_partiallycentralized_MPC_10 = 1.0e+04 * [2.713383194452476 2.700442585607222 2.696417323383893 2.686869527110244 2.682551331168068 2.670100215418421 2.665777484995863 2.664311680288444 2.657558688620042 2.664312277937123];
TFC_duplex_partiallycentralized_MPC_20 = 1.0e+04 * [2.712639958060689   2.699591247332050   2.693985517542243   2.689850285739326   2.678340228982182 2.673241209088883   2.664329134066781   2.660206585466814   2.659091883763914   2.660954143762251];
TFC_duplex_partiallycentralized_MPC_30 = 1.0e+04 * [2.715931378420213   2.699464562422975   2.692111718619767   2.688294094888666   2.680183987863412   2.671876118791451   2.662652184470312   2.673422903298604   2.658310474693437   2.666168092472544];
%TFC_duplex_partiallycentralized_MPC_40 = 1.0e+04 * [2.713383194452476 2.700442585607222 2.696417323383893 2.686869527110244 2.682551331168068 2.670100215418421 2.665777484995863 2.664311680288444 2.657558688620042 2.664312277937123];
%TFC_duplex_partiallycentralized_MPC_50 = 1.0e+04 * [2.713383194452476 2.700442585607222 2.696417323383893 2.686869527110244 2.682551331168068 2.670100215418421 2.665777484995863 2.664311680288444 2.657558688620042 2.664312277937123];

x = 1:10;

TFC_MPC_poly_parcent_5 = fit(x.',TFC_duplex_partiallycentralized_MPC_5.','poly2');
TFC_MPC_poly_parcent_10 = fit(x.',TFC_duplex_partiallycentralized_MPC_10.','poly2');
TFC_MPC_poly_parcent_20 = fit(x.',TFC_duplex_partiallycentralized_MPC_20.','poly2');
TFC_MPC_poly_parcent_30 = fit(x.',TFC_duplex_partiallycentralized_MPC_30.','poly2');


hold on
figure(1);
TFCnoMPC = (2.7647*1e4)*ones(1,length(TFC_duplex_partiallycentralized_MPC_30));
plot(TFCnoMPC, 'k', 'linewidth', 2);

plot(TFC_MPC_poly_parcent_5,'r');
plot(TFC_MPC_poly_parcent_10,'b');
plot(TFC_MPC_poly_parcent_20,'g');
%plot(TFC_MPC_poly_parcent_30,'k');

% title('TFC optimisations varying the number of vehicles','fontsize',14);
xlabel('Number of controlled CAVs','fontsize',15);
ylabel('Total fuel consumption [l/h]', 'fontsize',15);
xlim([1 10]);
% legend('TFC without MPC', 'TFC 1 hour optim. centralized', 'TFC with MPC optim. centralized','TFC 1 hour optim. decentralized', 'TFC with MPC optim. decentralized', 'TFC 1 hour optim. partially centralized', 'TFC with MPC optim. partially centralized');
legend('TFC without MPC', 'TFC MPC partially decent. radius 5','TFC MPC partially decent. radius 10','TFC MPC partially decent. radius 20', 'fontsize',15);
hold off

% % % % saveas(gca,'TFC_optimisations.png');
