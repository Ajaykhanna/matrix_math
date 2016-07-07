mat = A_mat_g;
myans = f_mat_g;

diagonal_check = true;
n = length(mat);
pivot_tol = 1.0e-5;

while(diagonal_check)
	diagonal_check = false;
	for i = 1:n
		if (mat(i,i) < pivot_tol) || (mat(i,i) > (-1)*pivot_tol)
			diagonal_check = true;
			mat(i,i) = 1.0;
			disp(i)
		end
	end
end