% -----Function description
function d = distance_b(C_matrix)
	% C_matrix
	% ---------c1;
	% ---------c2;
	% ---------c3;
	% .........ci;
	% ---------cm;

	[m, n] = size(C_matrix);
	d = 0;
	for i=1:m-1
		for j=i+1:m
			d_vectro = C_matrix(i, :) - C_matrix(j, :);
			d = d + norm(d_vectro, 2);
		end
	end