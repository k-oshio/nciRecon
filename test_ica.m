#!/usr/local/bin/octave

addpath("/Users/oshio/Math/ICA/fastica_2.1.octave");

# read input
x = csvread("IMG_x.csv");

# call fastICA
[s, A, W] = fastica(x, 'displayMode', 'off', 'verbose', 'off');
# output
csvwrite("IMG_s.csv", s);
csvwrite("ICA_a.csv", A);
