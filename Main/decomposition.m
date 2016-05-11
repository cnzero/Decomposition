function decomposition()
	clear,close all,clc
	channel = 1; % based on the single channel recognition.
	
	H = load('BasicInform.mat');
	%-Its head and tail have been cut.
	%-Three parts with the same movement have been connected.

	%----Database seems to be very necessary.----------


	%--Part one: Filtering with second order.
	%-rawdata
	rd = H.Data_EMG(channel, :);
	L = size(H.Data_EMG, 2);
	fd = [];
	for a=2:L-2
		fd = [fd, ...
			  rd(a+2)-rd(a+1)-rd(a)+rd(a-1)];
	end

	%- Compute the threshold value.
	C = 1.5;
	hold_value = C*sqrt(sum(fd.^2)/length(fd));

	ffd = [];
	for a=1:length(fd)
		if abs(fd(a))<hold_value
			ffd = [ffd, 0];
		else
			ffd = [ffd, fd(a)];
		end
	end

	L = floor(length(rd)/3);
	subplot(4,1,1)
	plot(rd(1:L))
	title('raw data of one channel')

	subplot(4,1,2)
	plot(fd(1:L))
	title('filted data with two orders')

	subplot(4,1,3)
	plot(ffd(1:L))
	title('filted data with threshold')

	% Peaklist(fd, hold_value);
	[Samples_x, Samples_y] = Peaklist(ffd, hold_value, 8);
    
%     plot some samples to check the results
    figure
    plot(Samples_y(10, :), 'r');
    hold on;
    plot(Samples_y(20,:), 'k');
    hold on;
    plot(Samples_y(30,:), 'b');
    hold on;
    plot(Samples_y(40,:), 'g');
    hold on;
    plot(Samples_y(50,:), 'y');

	% Mixture of Gaussian models. [GMM]
	% options = statset('Display', 'final');
	% for k=1:10
	% 	% k
	% 	try
	% 		obj = fitgmdist(Samples_y, k,'options', options);
	% 	catch exception
	% 		error = exception.message;
	% 		obj = fitgmdist(Samples_y, k,'options', options, 'Regularize', 0.1);
	% 	end
	% end













% function description
% to get the ultimate samples matrix
%    Peaks1
% 	 Peaks2
% 	 Peaks3
% 	 ......
% 	 Peaks_q  --->qxn

% Input:
%		the filted raw EMG data.
% Output:
% 		the ultimate samples matrix: x,y
function [MUAPx, MUAPy] = Peaklist(alist, hold_value, n)
	%-Mix the up and down peaks 
	[up_peak_index, down_peak_index] = PeaksCrossPoints(alist, hold_value);
	
	%--qx2
	%-mixing the upper and down peaks with the right order.
	peaks_index = sortrows([up_peak_index; ...
						   down_peak_index]);

	%-remove all zero-pair.
	% zero elements locates in the front part.
	% --1 method
	non_zero_index = 1;
	while(peaks_index(non_zero_index,1)==0)
		non_zero_index = non_zero_index + 1;
	end
	% slicing
	peaks_index = peaks_index(non_zero_index:size(peaks_index,1),:);

	MUAPx = []; % qx8
	MUAPy = [];	% qx8
	% from [strat_point end_point] to samples matrix.
	for i=1:size(peaks_index,1)
		Xs = Peaks_n(alist, peaks_index(i,1), peaks_index(i,2), n);
		MUAPx = [MUAPx;
				 Xs];
	end
	% overlapping cancle
    n_row = size(MUAPx, 1);
    MUAPx = Overlap(alist, MUAPx, n);
    new_n_row = size(MUAPx, 1);
    while ( n_row ~= (new_n_row+1) )
        n_row = size(MUAPx, 1);
        MUAPx = Overlap(alist, MUAPx, n);
        new_n_row = size(MUAPx, 1);
    end
	% ready to get MUAPy on position MUAPx
	for i=1:size(MUAPx, 1)
		y_row = [];
		for j=1:size(MUAPx, 2)
			y = alist(MUAPx(i,j));
			y_row = [y_row, y];
		end
		MUAPy = [MUAPy; ...
				 y_row];
	end

%------------------Function description
% find the upper peaks
% find the down peaks
%-Input: a sequent signal(t)
%		 hold_value, a positive threshold.
%-Output: 
% 	 x coordinates
%	 up_cross_p-->nx2[startp, endp], every row contains two cross point above hold_value
%	down_cross_p->nx2[startp, endp], every row contains two cross points under hold_value
function [up_cross_p, down_cross_p] = PeaksCrossPoints(alist, hold_value)
	L = length(alist);
	up_cross_p = zeros(L,2);
	down_cross_p = zeros(L,2);
	up_c = 1;
	dn_c = 1;

	for i=1:L-1
		d1 = alist(i) - hold_value;
		d2 = alist(i+1) - hold_value;

		d3 = alist(i) + hold_value;
		d4 = alist(i+1) + hold_value;

		if (d1<=0) && (d2>0) 		 %-A
			up_cross_p(up_c, 1) = i;

		elseif (d1>0) && (d2<=0)    %-B
			up_cross_p(up_c, 2) = i;
			up_c = up_c + 1;

		elseif (d3>0) && (d4<=0)     %-C
			down_cross_p(dn_c, 1) = i;

		elseif (d3<=0) && (d4>0)     %-D
			down_cross_p(dn_c, 2) = i;
			dn_c = dn_c + 1;
		end
	end

	%--Abstract the Non-zero part.
	up_cross_p = AbstractNonZero(up_cross_p);
	down_cross_p = AbstractNonZero(down_cross_p);


%-function description:
%-trim the set of all peaks list.
%-Input:
% 		the set of all peaks list
%-Output:
% 		the trimmed set of peaks list without zero index.
function newlist = AbstractNonZero(alist)
	%-alist
	% nx2
	if alist(1,1)==0
		startp = 2;
	else
		startp = 1;
	end
	count = size(alist, 1);
	if alist(count,1)==0 || alist(count,2)==0
		count = count -2;
	end

	newlist = alist(startp:count, :);

% function description:
% trim every Peaks list with length 8
% Input:
% 		all peaks lists
% Output:
% 		trimmed peaks x coordinate lists with length n.
function peak_array = Peaks_n(alist, startp, endp, n)
	max_value = alist(startp);
	max_index = startp;
	for i=startp+1:endp
		if alist(i)>max_value
			max_value = alist(i);
			max_index = i;
		end
	end
	peak_array = [];
	if mod(n, 2)==0
		% even_number
		sp = max_index - (n/2 - 1);
		ep = max_index + n/2;
	else
		% odd_number
		sp = max_index - floor(n/2);
		ep = max_index + floor(n/2);
	end
	peak_array = sp : ep;


% --------------Function description
% Input:
% 		Xs, index points of passing through threshold_value
% 				nx2
% 		n, 	the width of the peaks, generally 8 points
% Output:
% 		new_Xs, index points of passing through threshold_value
% 						 after trimming
%				mx2, m < n 					
function new_Xs = Overlap(alist, Xs, n)
	% part of matrix Xs
    % 9384        9385        9386        9387        9388        9389        9390        9391
    % 9430        9431        9432        9433        9434        9435        9436        9437
    % 9432        9433        9434        9435        9436        9437        9438        9439
    % 9438        9439        9440        9441        9442        9443        9444        9445
	new_Xs = [];
    i = 1;
	while ( i<=(size(Xs,1)-1) )
		% Overlap happends
		if (Xs(i, n)) > Xs(i+1, 1)
			new_x = Peaks_n(alist, Xs(i, 1), Xs(i+1, n), n);
            i = i+1;
        else
            new_x = Xs(i,:);
		end 
		new_Xs = [new_Xs; ...
				  new_x];
        i = i+1;
	end
