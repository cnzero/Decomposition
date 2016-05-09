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
	hold_value = 1*sqrt(sum(fd.^2)/length(fd));

	ffd = [];
	for a=1:length(fd)
		if abs(fd(a))<hold_value
			ffd = [ffd, 0];
		else
			ffd = [ffd, fd(a)];
		end
	end

	subplot(4,1,1)
	plot(rd)
	title('raw data')

	subplot(4,1,2)
	plot(fd)
	title('filted data with two orders')

	subplot(4,1,3)
	plot(ffd)
	title('filted data with threshold')

	% Peaklist(fd, hold_value);
	[Samples_x, Samples_y] = Peaklist(ffd, hold_value);

	plot_x = reshape(Samples_x', 1,[]);
	plot_y = reshape(Samples_y', 1,[]);
	subplot(4,1,4)
	% scatter(plot_x, plot_y, '.');
	title('plot all peaks ')

	% Mixture of Gaussian models. [GMM]
	options = statset('Display', 'final');
	for k=1:10
		k
		try
			obj = fitgmdist(Samples_y, k,'options', options)
		catch exception
			error = exception.message
			obj = fitgmdist(Samples_y, k,'options', options, 'Regularize', 0.1)
		end
	end

% function description
% to get hte ultimate samples matrix
%    Peaks1
% 	 Peaks2
% 	 Peaks3
% 	 ......
% 	 Peaks_q  --->qx8

% Input:
%		the filted raw EMG data.
% Output:
% 		the ultimate samples matrix: x,y
function [MUAPx, MUAPy] = Peaklist(alist, hold_value)
	%-Mix the up and down peaks 
	[up_peak_index, down_peak_index] = Peaks(alist, hold_value);
	
	%--qx2
	%-mixing the upper and down peaks with the right order.
	peaks_index = sortrows([up_peak_index; ...
						   down_peak_index]);

	%-remove all zero-pair.
	% zero elements locates in the front part.
	non_zero_index = 1;
	while(peaks_index(non_zero_index,1)==0)
		non_zero_index = non_zero_index + 1;
	end

	peaks_index = peaks_index(non_zero_index:size(peaks_index,1),:);

	MUAPx = []; % qx8
	MUAPy = [];	% qx8
	% from [strat_point end_point] to samples matrix.
	for i=1:size(peaks_index,1)
		Xs = Peaks_8(alist, peaks_index(i,1), peaks_index(i,2));
		Ys = [];

		for j=1:8
			Ys = [Ys, alist(Xs(j))];
		end
		MUAPx = [MUAPx;
				 Xs];
		MUAPy = [MUAPy;
				 Ys];
	end



%-function description
% find the upper peaks
% find the down peaks
%-Input: a sequent signal(t)
%		 hold_value, a positive threshold.
%-Output: 
% 	 x coordinates
%	 up_cross_p-->nx2[startp, endp], every row contains two cross point above hold_value
%	down_cross_p->nx2[startp, endp], every row contains two cross points under hold_value
function [up_cross_p, down_cross_p] = Peaks(alist, hold_value)
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
% 		trimmed peaks x coordinate lists with length 8.
function peak_array = Peaks_8(alist, startp, endp)
	max_value = alist(startp);
	max_index = startp;
	for i=startp+1:endp
		if alist(i)>max_value
			max_value = alist(i);
			max_index = i;
		end
	end
	peak_array = [max_index-3, ...
				  max_index-2, ...
				  max_index-1, ...
				  max_index, ...
				  max_index+1, ...
				  max_index+2, ...
				  max_index+3, ...
				  max_index+4];