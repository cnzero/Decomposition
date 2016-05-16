% function description

% Input:
% 		aList, raw sEMG data of one channel, 1xlength
% 		n_repeat, number of movement circle repeated.
% 		time_eachmove, how many second of each movement, 4s always
% 		n_cutsecond, how many second of tail or head slice of movement are cut, 0.4 for exampel

% Output:
% 		cell_movement, 
% 		cell_movement{1}, 1xa, for movement-1
% 		cell_movement{2}, 1xb, for movement-2
% 		. . . . .

function cell_movement = clipp_together(aList, n_movement, n_repeat, time_eachmove, n_cutsecond)
	cell_movement = cell(n_movement, 1);

	L = length(aList);

	% how many data points are acquired for one second.
	n_1second = floor(L/(n_movement - 1) /2/n_repeat/time_eachmove);
	% ~2000, the sampling frequency

	part_clipped = floor(n_1second * n_cutsecond);

	n_part = (n_movement - 1) * n_repeat * 2;

	for pt=0:n_part-2
		slice = aList(pt*time_eachmove*n_1second+part_clipped : ...
			          (pt+1)*time_eachmove*n_1second- part_clipped)
		if mod(pt, 2) == 0
			cell_movement{1} = [cell_movement{1}, slice];
		else
			seq = ceil(mod(pt, (n_movement - 1)*2)/2) + 1;
			cell_movement{seq} = [cell_movement{seq}, slice];
		end
	end