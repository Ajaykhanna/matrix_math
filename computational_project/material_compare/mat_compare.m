wm = fopen('mat_william.csv');
mat_william = textscan(wm,'%s');
mat_william = mat_william{1,:};
fclose(wm);

nick = fopen('mat_nick.csv');
mat_nick = textscan(nick,'%s');
mat_nick = mat_nick{1,:};
fclose(nick);

for i = 1:length(mat_william)
	if mat_william{i} ~= mat_nick{i}
		comparison(i) = {'DIFFERENT'};
	else
		comparison(i) = {' '};
	end
end