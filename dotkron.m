function X = dotkron(U,varargin)
% The Kroneck product of same rows for given matrices 
% More efficient than the first version DOTKRON
% E.g. U1 is I*J1, U2 is I*J2 and U3 is I*J3, then X is I*(J1 J2 J3)
%
% The main idea is referred to the code:
%  Laurent S (2020). Khatri-Rao product (https://www.mathworks.com/matlabcentral/fileexchange/28872-khatri-rao-product), 
%    MATLAB Central File Exchange. Retrieved June 18, 2020.

if ~iscell(U), U = [U varargin]; end
K = size(U{1},1);
if any(cellfun('size',U,1)-K)
    error('kr:RowMismatch', ...
          'Input matrices must have the same number of rows.');
end
I = size(U{1},2);
X = U{1};
for n = 2:length(U)
    J = size(U{n},2);
    A = reshape(U{n},[K 1 J]);
    X = reshape(bsxfun(@times,X,A),[K, I*J]);
    I = I*J;
end
X = reshape(X,[K size(X,2)]);
