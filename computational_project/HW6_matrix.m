clear
xs = define_xs;
N_g = length(xs.m87.f);


% infinte
% homogenous
% 8.7% MOX

total = eye(N_g);
for g = 1:N_g
	total(g,:) = total(g,:) * xs.m87.t(g);
end
scatter = xs.m87.s;
L = total - scatter;
L = L';

for g = 1:N_g
	for g_prime = 1:N_g
		P(g,g_prime) = xs.m87.x(g) * xs.m87.nu(g_prime) * xs.m87.f(g_prime);
	end
end

A = inv(L) * P;
k = eig(A);
[phi,~] = eig(A);

[keff, index] = max(k);
phi = abs(phi);
phi_fund = phi(:,2) / sum(phi(:,2));