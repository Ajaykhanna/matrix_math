close all


figure (1)
phi_mesh(1,:) = fund.phi(find(abs(mesh -  0.045) < 1e-5),:);
phi_mesh(2,:) = fund.phi(find(abs(mesh - 0.09) < 1e-5),:);
phi_mesh(3,:) = fund.phi(find(abs(mesh - 0.63) < 1e-5),:);
plot((1:N_g),phi_mesh(:,:));
legend('0.045','0.09','0.63')
xlabel('Group')
ylabel('phi')
print('test_b_1','-dpng','-r600');
close(1)

figure(2)
J_mesh(1,:) = fund.J(find(abs(mesh -  0.045) < 1e-5),:);
J_mesh(2,:) = fund.J(find(abs(mesh - 0.09) < 1e-5),:);
J_mesh(3,:) = fund.J(find(abs(mesh - 0.63) < 1e-5),:);

plot((1:N_g),J_mesh(:,:))
legend('0.045','0.09','0.63')
xlabel('Group')
ylabel('J')
print('test_b_2','-dpng','-r600');
close(2)

figure(3)
for i = 1:(length(mesh) - 1)
	average_x(i) = 0.5 * (mesh(i) + mesh(i + 1));
end
plot(average_x,fund.phibar((1:(length(average_x))),:));
legend('1','2','3','4','5','6','7')
xlabel('x [cm]')
ylabel('phi_bar','interpreter','none')
print('test_b_3','-dpng','-r600')
close(3)

figure(4)
for i = 1:size(fund.phibar,1)
	sum_phibar(i) = sum(fund.phibar(i,:));
end
plot(average_x,sum_phibar(:,(1:(length(average_x)))));
xlabel('x [cm]')
ylabel('sum(phi_bar)','interpreter','none')
print('test_b_4','-dpng','-r600')
close(4)

figure(5)
plot(mesh(1:(length(mesh) - 1)),fund.J((1:length(mesh) - 1),:))
legend('1','2','3','4','5','6','7')
xlabel('x [cm]')
ylabel('J')
print('test_b_5','-dpng','-r600')
close(5)
