mat = A_mat_g;
myans = f_mat_g;

diagonal_check = true;
n = length(mat);
pivot_tol = 1.0e-5;
while (diagonal_check)
	diagonal_check = false;
	for i = 1:n
		if ((mat(i,i) < pivot_tol) && (mat(i,i) > (-1)*pivot_tol))
			fprintf('pivot  %i \n',i);
			diagonal_check = true;
			store_row = mat(i,:);
			store_myans = myans(i);
			
			% find a swap candidate
			it = 0;
			j = i;
			while (it <= n)
				it = it + 1;
				j = j - 1;
				if (j < 1)
					j = ;
				end
				if ((mat(j,i) > pivot_tol) || (mat(j,i) < (-1)*pivot_tol))
					% if ((mat(i,j) > pivot_tol) || (mat(i,j) < (-1)*pivot_tol))
						mat(i,:) = mat(j,:);
						myans(i) = myans(j);
						mat(j,:) = store_row;
						myans(j) = store_myans;
						fprintf('swap  %i \n',j);
						break
					% end
				end
			end
		end
	end
end
