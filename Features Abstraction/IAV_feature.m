%[IAV-Integrated Absolute Value] feature
function feature = IAV_feature(data_TimeWindow)
	feature = sum(abs(data_TimeWindow));
