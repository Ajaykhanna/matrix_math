clear
%
% USER INPUT
%

% test a, b, or c
test_char = 'c';

% import test geometry and cross sections
run(sprintf('test_%s.m',test_char));
xs = define_xs; % external function
N_g = length(xs.uo2.t);

%
% INPUT PROCESSING
%

% check for proper number of materials
for i = 1:length(pin)
	if length(pin(i).x) ~= length(pin(i).mat)
		error('pin.x doesnt match pin.mat')
	end
end
% make big list of x-coordinates and materials for the entire problem
for i = 1:length(pin_map)
	if i == 1
		x = pin(pin_map(i)).x;
		mat = pin(pin_map(i)).mat;
	else
		x = [x (x(end) + pin(pin_map(i)).x)];
		mat = [mat pin(pin_map(i)).mat];		
	end
end
capital_H = x(end);
mesh = dx:dx:x(end);
% determine what material/cross sections at each mesh point
mesh_tol = 1e-3;
for i = 1:length(mesh)
	region(i).x = mesh(i);
	for j = 1:length(x)
		if (region(i).x - x(j)) <= mesh_tol
			region(i).mat = mat{j};
			region(i).xs = eval(sprintf('xs.%s',region(i).mat));
			break
		end
	end
end

% check sum of chi for fissionable materials
chi_tol = 1e-5;
chi_check = false;
if abs(sum(xs.uo2.x) - 1) >= chi_tol
	chi_check = true;
end
if abs(sum(xs.m43.x) - 1) >= chi_tol
	chi_check = true;
end
if abs(sum(xs.m70.x) - 1) >= chi_tol
	chi_check = true;
end
if abs(sum(xs.m87.x) - 1) >= chi_tol
	chi_check = true;
end
if chi_check
	error('chi check failed.');
end

%
% CALCULATIONS
%

% build A matrix of geometry/material properties
A = zeros((3 * N + 2),(3 * N + 2),N_g);
for g = 1:N_g
	for i = 1:N
		A((3 * i - 2),(3 * i - 1),g) = -1;
		A((3 * i - 2),(3 * i),g) = region(i).xs.r(g) * dx;
		A((3 * i - 2),(3 * i + 2),g) = 1;
		
		A((3 * i - 1),(3 * i),g) = (-1) * region(i).xs.D(g);
		A((3 * i - 1),(3 * i + 1),g) = region(i).xs.D(g);
		A((3 * i - 1),(3 * i + 2),g) = 0.5 * dx;
		
		A((3 * i),(3 * i - 2),g) = (-1) * region(i).xs.D(g);
		A((3 * i),(3 * i - 1),g) = 0.5 * dx;
		A((3 * i),(3 * i),g) = region(i).xs.D(g);
	end
end
% reflective BC in A matrix
A((3 * N + 1),1,:) = -1;
A((3 * N + 1),3,:) = 1;
A((3 * N + 2),(3 * N + 1),:) = -1;
A((3 * N + 2),(3 * N),:) = 1;


% initial guesses
s = 1;
k = 1;
phibar(:,:,s) = (ones(N,N_g,1)) * (1/ (N * N_g));
phibarerror = 1;
f = zeros((3 * N + 2),N_g);

% convergence tolerances
epsilon_k = 1e-5;
epsilon_phi = 1e-5;

while phibarerror > epsilon_phi || kerror > epsilon_k 
	s = s + 1;
	% check for number of iterations
	if s == 101;
		error('failure to converge after 100 iterations')
	end
	% calculate Q values
	for g = 1:N_g
		for i = 1:N
			fission_term = 0;
			for g_prime = 1:N_g
				fission_term = fission_term + region(i).xs.nu(g_prime) * region(i).xs.f(g_prime) * phibar(i,g_prime,(s - 1));
			end
			Q_f(i,g,(s - 1)) = (region(i).xs.x(g) / k(s - 1)) * fission_term;

			upscatter_term = 0;
			for g_prime = (g + 1):N_g
				upscatter_term = upscatter_term + region(i).xs.s(g_prime,g) * phibar(i,g_prime,(s - 1));
			end
			Q_up(i,g,(s - 1)) = upscatter_term;

			downscatter_term = 0;
			for g_prime = 1:(g - 1)
				downscatter_term = downscatter_term + region(i).xs.s(g_prime,g) * phibar(i,g_prime,s);
			end
			Q_down(i,g,s) = downscatter_term;
		
			Q(i,g,s) = Q_down(i,g,s) + Q_up(i,g,(s - 1)) + Q_f(i,g,(s - 1));
			% turn Q's into f vector
			f((3 * i - 2),g) = Q(i,g,s) * dx;
		end
		
		% perform matrix math to sovle for y vector
		y(:,g) = A(:,:,g) \ f(:,g);
		
		% parse y vector for required functions
		for i = 1:(N + 1)
			if i <= N
				phi(i,g,s) = y((3 * i - 2),g);
				J(i,g,s) = y((3 * i - 1),g);
				phibar(i,g,s) = y((3 * i),g);
			else
				phi(i,g,s) = y((3 * i - 2),g);
				J(i,g,s) = y((3 * i - 2),g);
			end
		end
		
	end
	
	% update the guess for k
	numerator = 0;
	denominator = 0;
	for g_prime = 1:N_g
		for i = 1:N
			numerator   = numerator   + region(i).xs.nu(g_prime) * region(i).xs.f(g_prime) * phibar(i,g_prime,s) * dx;
			denominator = denominator + region(i).xs.nu(g_prime) * region(i).xs.f(g_prime) * phibar(i,g_prime,(s - 1)) * dx;
		end
	end
	
	% make sure k is positive
	k(s) = k(s - 1) * (numerator / denominator);
	if k(s) < 0
		error('negative k value')
	end
	
	% normalized calculated funcitons
	sum_phibar = 0;
	sum_phi    = 0;
	for g = 1:N_g	
		for i = 1:N
			sum_phibar = sum_phibar + phibar(i,g,s) * dx;
			sum_phi    = sum_phi    + phi(i,g,s) * dx;
		end
	end			
	phibar(:,:,s) = (phibar(:,:,s) * capital_H) / (sum_phibar);
	phi(:,:,s)    = (phi(:,:,s)    * capital_H) / (sum_phi);
	
	% calculate iteration error for the next iteration
	phibarerror = (norm(phibar(:,:,s) - phibar(:,:,(s - 1)),Inf)) / norm(phibar(:,:,s),Inf);
	kerror = abs(1 - (k(s - 1) / k(s)));
end

% save the fundamental mode to a variable
fund.phibar = phi(:,:,size(phibar,3));
fund.phi = phi(:,:,size(phi,3));
fund.J = J(:,:,size(J,3));
fund.k = k(end);

%
% OUTPUT
%

% run script to generate plots
run(sprintf('data_%s.m',test_char));
% output numbers to file
fname = fopen(sprintf('test_%s.out',test_char),'w');
fprintf(fname,'Number of iterations: %i \n',s);
fprintf(fname,'keff: %.6f',fund.k);