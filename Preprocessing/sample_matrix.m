% functional description
% Input:
% 		aList, mxn, 
% 				the rawdata of one movement--with clipped together
% 				m, number of channels
% 				n, length of each movement
% 		LW, the length of time window, 128 for example
% 		LI, the length of increase window, 64 for example
% 		features_name, {'SSC', 'ZC', 'IAV'}, for example.
% 		
% Output:
% 		Sample_Matrix
% 				n_window x (n_ch x n_f)
% 		p1:	ch1f1, ch1f2, ch1f3,...,ch2f1, ch2f2, ch2f3...
% 		p2: ch1f1, ch1f2, ch1f3,...,ch2f1, ch2f2, ch2f3...
% 		. . . . . .

function Matrix = sample_matrix(aList, LW, LI, features_name)
	Matrix = [];
	n_channel = size(aList, 1);
	n_window = floor(1 + ( size(aList, 2) - LW )/LI);

	for nw=0:n_window-1
		row = [];
		for ch=1:n_channel
			feature = [];
			data = aList(ch, nw*LI+1 : nw*LI+LW);
			for fn=1:length(features_name)
				func_handle = str2func(features_name{fn});
				feature = [feature, func_handle(data)];
			end
			row = [row, feature];
		end
		Matrix = [Matrix; ...
					  row];
	end

	

% Realization of useful feature extraction.
% -------------Attention-------
% some functions need more than one input argument

% features list
% ARC
% Ceps
% IAV
% LogD
% MAV
% MDF
% MNF
% PSD, no realization
% RMS
% SSC
% VAR
% WA
% WL
% ZC
% MAX
% MED
% Non_MED
% SemiEny1
% SemiEny2

% ARC, Auto Regression Model coefficients.
% the first coefficient of ARC with order four.
function feature = ARC(x)  
	order = 4;              
	cur_xlpc = real(lpc(x,order)');
	feature = -cur_xlpc(order+1,:);

% Cepstrum coefficients, Ceps
% Equation(1-11) in Xiong's Ph.D article.
function feature = Ceps(data_TimeWindow, order)
	arc = getar_lpcfeat(data_TimeWindow, order);
	%size: orderx1
	
	feature(1, 1) = -arc(1);
	for l=2:order
		total = 0;
		for j=1:l-1
			total = total + (1 - j/l)*arc(j)*feature(1, l-j);
			feature(1, l) = -arc(l) - total;
		end
	end
function feat = getar_lpcfeat(x,order)                %the order of parameter in feat：[a1; a2; a3; a4]
	cur_xlpc = real(lpc(x,order)');
	feat = -cur_xlpc(2:(order+1),:);
%size: orderx1

%[IAV-Integrated Absolute Value] feature
function feature = IAV(data_TimeWindow)
	feature = sum(abs(data_TimeWindow));

function feature = LogD(data_TimeWindow)
	feature = exp(mean(log(abs(data_TimeWindow))));

% [MAV-Mean Absolute Value] feature
function feature = MAV(data_TimeWindow)
	feature = mean(abs(data_TimeWindow));

% [MDF-Median Frequency] feature
function feature = MDF(data_TimeWindow)
	Window_L = length(data_TimeWindow);
	power_spectral = abs(fft(data_TimeWindow, Window_L).^2)/Window_L;
	L = length(power_spectral);
	Fs = 2000; % Sampling frequency
	f = Fs/2*linspace(0, 1, L/2)';
	feature = medianfrequency(f, power_spectral(1:L/2));
% bisection search for increased speed.
function mnf = medianfrequency(f,p)
	N = length(f);
	low = 1;
	high = N;
	mid = ceil((low+high)/2);
	while ~( sum(p(1:mid)) >= sum(p((mid+1):N)) && sum(p(1:(mid-1))) < sum(p(mid:N)) )
	    if sum(p(1:mid)) < sum(p((mid+1):N))
	        low = mid;
	    else
	        high = mid;
	    end
		% fixed the original bug here.
        if p(1) == max(p)
            mid = floor((low+high)/2);
        else
            mid = ceil((low+high)/2);
        end
    end
	mnf = f(mid);

% Mean Frequency, MNF
function feature = MNF(data_TimeWindow)
	Window_L = length(data_TimeWindow);
	power_spectral = abs(fft(data_TimeWindow, Window_L).^2)/Window_L;
	L = length(power_spectral);
	Fs = 2000; % Sampling frequency
	f = Fs/2*linspace(0, 1, L/2)';
	feature = meanfrequency(f, power_spectral(1:L/2));
function mnf = meanfrequency(f, p)
	mnf = sum(p.*f)/sum(p);

% Power Spectral Density, PSD
function feature = PSD(data)
	% 918

% [RMS-Root Mean Square]
function feature = RMS(data_TimeWindow)
	feature = sqrt(mean(data_TimeWindow.^2));

%[SSC-Slope Sign Change] feature
function feature = SSC(data_TimeWindow, DeadZone)
	data_size = length(data_TimeWindow);
	feature = 0;

	if data_size == 0
		feature = 0;
	else
		for j=3:data_size
			difference1 = data_TimeWindow(j-1) - data_TimeWindow(j-2);
			difference2 = data_TimeWindow(j-1) - data_TimeWindow(j);
			Sign = difference1 * difference2;
			if Sign > 0
				if abs(difference1)>DeadZone || abs(difference2)>DeadZone
					feature = feature + 1;
				end
			end
		end
		feature = feature/data_size;
	end


% Variance - VAR
function feature = VAR(data_TimeWindow)
	feature = var(data_TimeWindow);

%-[WA-Willison Amplitude] feature
function feature = WA(data_TimeWindow, Threshold)
	data_size = length(data_TimeWindow);
	feature = 0;
	if data_size == 0
		feature = 0;
	else
		for i=2:data_size
			difference = data_TimeWindow(i) - data_TimeWindow(i-1);
			if abs(difference)>Threshold
				feature = feature + 1;
			end
		end
		feature = feature/data_size;
	end


%[WL-Waveform Length] feature
function feature = WL(data_TimeWindow)
	feature = sum(abs(diff(data_TimeWindow)))/length(data_TimeWindow);

%[ZC-zero crossing] feature
function feature = ZC(data_TimeWindow, DeadZone)
	data_size = length(data_TimeWindow);
	feature = 0;

	if data_size == 0
		feature = 0;
	else
		for i=2:data_size
			difference = data_TimeWindow(i) - data_TimeWindow(i-1);
			multy      = data_TimeWindow(i) * data_TimeWindow(i-1);
			if abs(difference)>DeadZone && multy<0
				feature = feature + 1;
			end
		end
		feature = feature/data_size;
	end


% [MAX] maximium value feature
function feature = MAX(data)
	feature = max(abs(data));

% [MED] median value feature
function feature = MED(data)
	feature = median(data); 

% [Non_MED] 
% ------Attention----
% NaN value return
function feature = Non_MED(data)
	feature = median(data(data~=0));

% SemiEny1
function feature = SemiEny1(data)
	N = length(data);
	feature = sum(data(floor(N/2)+1:end).^2);

% SemiEny2
function feature = SemiEny2(data)
	N = length(data);
	feature = sum(data(1:floor(N/2)).^2);
