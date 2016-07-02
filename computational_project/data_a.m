figure (1)
phi_mesh(1,:) = fund.phi(find(abs(mesh -  0.54) < 1e-5),:);
plot((1:N_g),phi_mesh(:,:));
xlabel('Group')
ylabel('phi')
print('test_a_1','-dpng','-r600');
close(1)