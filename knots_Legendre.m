% Gauss-Legendre quadrature rule on [x_a,x_b]
% n: # of quadrature points (integer)
% x_a, x_b: inf and sup of domain (x_a<x_b)
% whichrho: 'prob' means that weight function is 1/(x_b-x_a),
% 'nonprob' means that weight function is 1. 
% Default whichrho is 'prob'
% x: quadrature points (1 by n vector)
% w: quadrature weights (1 by n vector)
% x_a=-1, x_b=1, whichrho='nonprob' gives standard Gauss-Legendre quadrature rule
% Matched with lev2knots_lin

function [x,w]=knots_Legendre(n,x_a,x_b,whichrho)

if nargin==3
    whichrho='prob';
end

if n==1

    x=(x_a+x_b)/2;
    wt=1;

else
    
    % calculates the values of the recursive relation
    [a,b]=coeflege(n);
    
    % builds the matrix
    JacM=diag(a)+diag(sqrt(b(2:n)),1)+diag(sqrt(b(2:n)),-1);
    
    % calculates points and weights from eigenvalues / eigenvectors of JacM
    [W,X]=eig(JacM);
    x=diag(X)';
    wt=W(1,:).^2;
    [x,ind]=sort(x);  %#ok<TRSRT>
    wt=wt(ind);
    
    % modifies points according to the distribution and its interval x_a, x_b
    x = (x_b-x_a)/2*x + (x_a+x_b)/2;
    
end

% finally, fix weights

switch whichrho
    
    case 'nonprob'
        w=(x_b-x_a)*wt;

    case 'prob'
        w=wt;

    otherwise
    error('SparseGKit:WrongInput','4th input not recognized')
    
end



%----------------------------------------------------------------------
function [a, b] = coeflege(n)

if (n <= 1), disp(' n must be > 1 '); 
    return; 
end

a=zeros(1,n);
b=zeros(1,n);

b(1)=2;

k=2:n;
b(k)=1./(4-1./(k-1).^2); 