function out = upperTri(X)

% Either:
%    - for matrix X, extract and vectorise the upper triangular part.
%    - for vector X with 1/2*n(n+1) elements (n an integer) then form an
%        upper triangular matrix
%
%  i.e. note that upperTri(upperTri(X)) is X with the lower-triangular part
%  replaced with zeros

siz = size(X);

if length(siz) == 2 && siz(1) == siz(2)
    out = X(triu(ones(siz(1)))==1);
elseif length(siz) == 2 && (siz(1) == 1 || siz(2) == 1)
    p = sqrt(1/4+2*numel(X))-1/2;
    if p == floor(p)
        out = zeros(p);
        for k = 1:p
            out(1:k,k) = X(1/2*k*(k-1)+1:1/2*k*(k+1));
        end
    else
        error('Vector input X must have a number of elements equal to a triangular number');
    end
else
    error('Input X not appropriately size - check function spec');
end