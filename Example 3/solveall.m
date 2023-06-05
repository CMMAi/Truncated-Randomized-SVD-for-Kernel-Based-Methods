function [errintmat, errevalmat]...
    =solveall(Aint, normAint, Pint, Aeval, Peval, findices, methindices)
global raTa
% applies all solvers to all examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results are matrices, rows = examples, columns = solvers
nm=length(methindices); % the actual list of solvers requested
ralist=zeros(nm,1);
tolElist=zeros(nm,1);
% we check first if some SVD is needed
for k=1:length(methindices)
    kmeth=methindices(k);
    switch(kmeth)
        case 1
            [U, S, V]=svd(Aint);SD=diag(S);
            tol = max(size(Aint)) * eps(max(SD));
            ralist(k,1)=sum(SD > tol);
%             tolElist(k,1)  = RSeps(S(1,1));
        case 2
            [Ur, Sr, Vr]=rsvd0(Aint,raTa);SrD=diag(Sr);
%             ralist(k,1)=RSrank(Sr);
            tol = max(size(Aint)) * eps(max(SrD));
            ralist(k,1)=sum( SrD > tol);
            %             tolElist(k,1)  = RSeps(Sr(1,1));
    end
end

for i=1:length(findices) % run over all functions in list
    fin=findices(i); % get the function number
    beval=myfct(Peval,fin); % fct on evaluation points
    rng(1)
    bint = myfct(Pint, fin)+ 0.000*randn(size(Pint,1),1); % on interpolation points
    normbint=norm(bint,Inf); % norm
    % now run over solvers, with somewhat different interfaces.
    % Add other solvers here. 
    for k=1:length(methindices)
        kmeth=methindices(k);
        switch(kmeth)
            case 1   
                [errint, erreval]...
                 =tSVD(Aint, bint, normAint, ralist(k), U, S, V, normbint, Aeval, beval);
            case 2   
                [errint, erreval]...
                 =trSVD(Aint, bint, normAint, ralist(k), Ur, Sr, Vr, normbint, Aeval, beval);
            otherwise
                error('unimplemented method')
        end
        errintmat(i,k)=errint;
        errevalmat(i,k)=erreval;

    end
end

function [errint, erreval, coeff]...
          =tSVD(Aint, bint, normAint, ra, U, S, V, normbint, Aeval, beval)
% the solver for SVD truncated at rank
c=U'*bint;
y=S\c;
coeff=V(:,1:ra)*y(1:ra,1);
errint= norm(bint - Aint *coeff, Inf);
erreval=norm(beval-Aeval*coeff, Inf);

function  [errint, erreval, coeff]...
           =trSVD(Aint, bint, normAint, ra, U, S, V, normbint, Aeval, beval)
% the solver for SVD truncated at rank
c=U'*bint;
y=S\c;
coeff=V(:,1:ra)*y(1:ra,1);
errint= norm(bint - Aint *coeff, Inf);
erreval=norm(beval-Aeval*coeff, Inf);
