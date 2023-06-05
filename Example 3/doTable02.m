function doTable02()
clear all;
close all;
% change RBFs: For Gaussian: type:'g' with par=any values
%              FOR MQ:  type='mq' with par= 0.5
%              FOR IMQ: type='mq' with par= -0.5
%              For W2:  type='w' with par= 2 **see below note.
%              For Matern with \nu=5:  type='ms' with par= 5
type='g'; 
par=1;
% n: number of tentative interpolation points
n=1000;
% nt:  number of  evaluation points
nt=900;
% findices: Test function indices: "1" for F1 and "2" for F2
findices=[1,2];
% method: method indices: "1" for tSVD and "2" for trSVD
method = [1,2];
errList(type, par, n, nt, findices, method)

%** note on Wendland RBF: When using Wendland, change "scaleList" in 
%"Ex3_errList" function to "scaleList = RBFscalestart:0.1:20;"
%since the Wendland RBF needs higher scale values to converge.
% other RBFs are fine with default "scaleList = RBFscalestart:0.01:1;"

%%% main code and instruction can be found on 
% "https://num.math.uni-goettingen.de/schaback/research/papers/MPfKBM.pdf"

