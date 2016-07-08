mat = A_mat_g;
myans = f_mat_g;

diagonal_check = true;
n = length(mat);
pivot_tol = 1.0e-5;
non_zero = zeros(n,1);

for i = 1:n
	% non_zero(i) = 0;
	for j = 1:n
		if (mat(i,j) > pivot_tol) || (mat(i,j) < (-1)*pivot_tol)
			non_zero(i) = non_zero(i) + 1;
		end
	end
end

z = 0;
diagonal_check = true;
while (diagonal_check)
	done = false(n);
	diagonal_check = false;
	z = z + 1;
	for x = min(non_zero):max(non_zero)
		for i = 1:n
			if (non_zero(i) == x) && ((mat(i,i) < pivot_tol) && (mat(i,i) > (-1)*pivot_tol))
				diagonal_check = true;
				for j = 1:n
					if (mat(i,j) > pivot_tol) || (mat(i,j) < (-1)*pivot_tol)
						if (not(done(j)))
							mat([i j],:) = mat([j i],:);
							done(j) = true;
						end
					end
				end
			end
		end
	end
end

disp(z)
for i = 1:n
	if (mat(i,i) < pivot_tol) && (mat(i,i) > (-1)*pivot_tol)
		disp(i)
		error('you were wrong')
	end
end