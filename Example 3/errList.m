function Ex3_errList(type, par, n, nt, findices, method)
% Checking beyond the Condition Limit but in full-rank setting
% nt:  number of  evaluation points
% n: number of tentative interpolation points
% findices: Test function indices: 1 for F1 and 2 for F2
% method: method indices: 1: for tSVD and 2 for trSVD
close all
global RBFscale
global RBFpar %=par 
global RBFtype %=type
global raTa
RBFtype=type
RBFpar=par
RBFscale=1
raTa = n/2; % rank target for trSVD algorithm
nc=1*n; % number of  center points
% findices all functions we want, specify a list if not all are wanted
a=0;b=1;  %Domain [a b]^2 
[Pint, Pcntr, Peval, Pex, Pey]=getPoints01(n,nc,nt,a,b);
RBFscalestart=0.01;%0.1/sqrt(n) % 0.5;my starting scale for [0,1]^2
nPint=length(Pint(:,1))
nPeval=length(Peval(:,1)) % evaluation points
% lists of output matrices
errintmatlist=[]; % interpolation error
errevalmatlist=[]; % evaluation error 
methindices=[method];
scaleList = RBFscalestart:0.01:1;
% loop over scales
for RBFscale=scaleList
    % get matrices
    Aint=kermat(Pint,Pcntr);
    Aint=Aint+(rand(n,nc)*2-1)*0.000; % add noise % @@
    Aeval=kermat(Peval,Pcntr); 
    % [co, ra, sra, mycond, SVDrepro, SVDreprora,U,S,V]=MatrixExam(Aint,'plot');
    normAint =norm(Aint,Inf); % we need norms later for the evaluation condition 
    % solve for all functions and various methods
    [errintmat, errevalmat]...
        =solveall(Aint, normAint, Pint, Aeval, Peval, findices, methindices);
    % rank check, save scale if rank still OK

    % now add results. Each mat is a matrix functions times methods. 
    errintmatlist=[errintmatlist errintmat];
    errevalmatlist=[errevalmatlist errevalmat];
    RBFscale % just for seeing if the program halts
end

error_tSVD=[];
error_trSVD=[];
disp('F,  &c(trsvd)  &error(trsvd)  &c(tsvd)    &error(tsvd)')
for fid=1:length(findices)
    fin=findices(fid);
    
    error_tSVD(fid,:)  = errevalmatlist(fid,1:2:end);
    [min_tsvd,idC_tsvd] = min(error_tSVD(fid,:));
    
    error_trSVD(fid,:) = errevalmatlist(fid,2:2:end);
    [min_trsvd,idC_trsvd] = min(error_trSVD(fid,:));
    
    disp(sprintf('%d    &%3.3f      &%2.2e       &%3.3f       &%2.2e   \\\\',fin,  scaleList(idC_tsvd), min_tsvd, scaleList(idC_trsvd), min_trsvd))
    
end
end



