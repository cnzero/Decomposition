y = decomposition();
close all,clc;

for c=2:10
	c
	[idx, C, sumd] = kmeans(y,c);
	distance_within = sum(sumd);
	distance_between = distance_b(C);
	LDA_value = distance_between/distance_within
end
