function I = binsrchrightn(X,F) %#codegen
%BINSRCHRIGHTN  Binary search for rightmost insertion index, many times.
%
%
% I = binsrchrightn(X,F) 
% finds the indices I(n) for which
%   F(I(n),n) <= X(n) < F(I(n)+1,n).
% given a vector X and a matrix F with sorted columns. All NaN's in F
% behave as if they were +Inf. When X(n) < F(1,n), then I(n) = 1. When X(n)
% > F(end,n), then I(n) = size(F,1). If X is a scalar, it is used in place
% of X(n) above. If F is a sorted column vector, F(j) replaces F(j,n)
% above.
%
% 
% This function is very similar to MATLAB's discretize, except 
% (a) that binsrchrightn allows F to be multi-dimensional, whereas
%     discretize requires its 'edges' input to be a vector, and
% (b) binsrchrightn returns 1 or size(F,1) where X is outside the range of
%     F, whereas bdiscretize returns NaN here.
%
%
% --- Input:
% X [vector of length N]: evaluation sites
% F [M,N]: data sites
%
%
% --- Output:
% I [vector or length N]: indices
%
%
% --- Code generation:
% A flexible use of codegen is as follows, given some values for nx, ny,
% and nz:
%     mexconfig = coder.config('mex');
%     mexconfig.ExtrinsicCalls = false;
%     mexconfig.ResponsivenessChecks = false;
%     mexconfig.IntegrityChecks = false;
%     type1 = coder.typeof(0, [nx, ny], [true true]);
%     type2 = coder.typeof(0, [nz, nx, ny], [false, true true]);
%     codegen('binsrchrightn', '-o', 'binsrchrightn_mex', '-args', {type1, type2}, '-config', mexconfig);
% Note: in MATLAB 2015b and above (with the JIT compiler), the uncompiled 
%       version is not much slower than the codegen version.
%
%
% Author    : Geoff Stanley
% Email     : geoffstanley@gmail.com
% Version   : 1.0
% History   : 17/01/2019 - initial release


szF = size(F);
FN = prod(szF(2:end));
K = szF(1);

% Decide which inputs to loop through and which to stick with 1st column
if numel(X) == FN
    N = FN;
    I = squeeze(zeros(szF(2:end)));
    nK = 0;
    for n = 1:N % Loop over columns, each of which is a binary search problem
        x = X(n);
        i = 1; % Result will be >= 1 always. May be useful to change i=2, e.g. for linear interpolation.
        k = K; % Result will be <= K always
        while (i < k)
            j = ceil((i + k) / 2);
            if F(j+nK) <= x
                i = j;
            else
                k = j - 1;
            end
        end
        I(n) = i;
        nK = nK + K;
    end
elseif FN == 1
        N = numel(X);
        I = squeeze(zeros(size(X)));
        for n = 1:N % Loop over columns, each of which is a binary search problem
            x = X(n);
            i = 1; % Result will be >= 1 always. May be useful to change i=2, e.g. for linear interpolation.
            k = K; % Result will be <= K always
            while (i < k)
                j = ceil((i + k) / 2);
                if F(j) <= x
                    i = j;
                else
                    k = j - 1;
                end
            end
            I(n) = i;
        end
else % numel(X) == 1
    N = FN;
    I = squeeze(zeros([1, szF(2:end)]));
    nK = 0;
    for n = 1:N % Loop over columns, each of which is a binary search problem
        i = 1; % Result will be >= 1 always. May be useful to change i=2, e.g. for linear interpolation.
        k = K; % Result will be <= K always
        while (i < k)
            j = ceil((i + k) / 2);
            if F(j+nK) <= X
                i = j;
            else
                k = j - 1;
            end
        end
        I(n) = i;
        nK = nK + K;
    end
end
end