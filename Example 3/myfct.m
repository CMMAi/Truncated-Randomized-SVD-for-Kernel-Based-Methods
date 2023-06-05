function fval=myfct(Points, ind)
% the list of unctions for examples
x=Points(:,1);
y=Points(:,2);
switch ind
    case 1
        fval=(1/9)*(tanh(9*y-9*x)+1);
    case 2
        fval=0.0025./((x-1.01).^2+(y-1.01).^2);
    otherwise
        error('Case not implemented')
end