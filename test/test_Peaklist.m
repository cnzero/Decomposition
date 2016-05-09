function test_Peaklist()
	%-preprocessing
	close all,clear,clc
	x = 0:0.01:6*pi;
	y = -cos(x);
	plot(x,y,'-.')
	[up, down] = Peaklist(y,0.8)

	for i=1:length(up)
		hold on;
		plot(x(up(i,1):up(i,2)),y(up(i,1):up(i,2)), 'r');
	end

	for i=1:length(down)
		hold on;
		plot(x(down(i,1):down(i,2)),y(down(i,1):down(i,2)), 'k');
	end




%-function description
%-Input: a sequent signal(t)
%		 hold_value, a positive threshold.
%-Output: 
%	 up_cross_p-->nx2, every row contains two cross point above hold_value
%	down_cross_p->nx2, every row contains two cross points under hold_value
function [up_cross_p, down_cross_p] = Peaklist(alist, hold_value)
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
		% save('up_cross_p')
	end

	%--Abstract the Non-zero part.
	up_cross_p = AbstractNonZero(up_cross_p, up_c);
	down_cross_p = AbstractNonZero(down_cross_p, dn_c);


function newlist = AbstractNonZero(alist, count)
	%-alist
	% nx2
	if alist(1,1)==0
		startp = 2;
	else
		startp = 1;
	end
	if alist(count,1)==0 || alist(count,2)==0
		count = count -1;
	end

	newlist = alist(startp:count, :);
