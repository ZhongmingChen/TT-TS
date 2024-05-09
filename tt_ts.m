function [G, err, rel, T] = tt_ts(Y, R, J1, sigma, varargin)
% TT_TS     Implementation of TT-TS algorithm. TT-TS
%           utilizes TensorSketch to compute the Tensor train decomposition
%           of a tensor.
%
%           This function requires Tensor Toolbox [1] version 2.6. 
% 
%   [G, err, rel, T] = tt_ts(Y, R, J1, sigma, varargin) returns an approximate rank-R TT
%   decomposition of Y in the form of a core tensor G. Y can be a Matlab double array or one of the following types from 
%   Tensor Toolbox: ktensor, tensor, ttensor, or sptensor. R is a vector 
%   containing the target dimension, and J1 is the  sketch 
%   dimensions used. G is a Tensor Toolbox tensor, err is the error of 
%   core, rel is the error of tensor, and T is the time of each sweep. 
%
% REFERENCES:
%
%   [1] B. W. Bader, T. G. Kolda and others. MATLAB Tensor Toolbox 
%       Version 2.6, Available online, February 2015. 
%       URL: http://www.sandia.gov/~tgkolda/TensorToolbox/.
%
%   %[2] O. A. Malik, S. Becker. Low-Rank Tucker Decomposition of Large 
%       Tensors Using TensorSketch. Advances in Neural Information 
%       Processing Systems (NeurIPS), 2018.

%% Include relevant files

addpath(genpath('help_functions'));

%% Handle optional inputs

params = inputParser;
addParameter(params, 'tol',1e-12, @isscalar);
addParameter(params, 'maxiters', 100, @(x) isscalar(x) & x > 0);
addParameter(params, 'verbose', false, @isscalar);
parse(params, varargin{:});
tol = params.Results.tol;
maxiters = params.Results.maxiters;
verbose = params.Results.verbose;

%% Set extflag
% extflag will be true if we are using an external function to compute the sketches, otherwise it will be false

    if iscell(Y)
        extflag = true;
    else
        extflag = false;
    end

%% Convert input tensor to double array if necessary. Set sflag (sparsity flag)

    sflag = false;
    if ismember(class(Y), {'ktensor', 'tensor', 'ttensor'})
        Y = double(Y);
    elseif ismember(class(Y), {'sptensor'})
        sflag = true;
    elseif ~ismember(class(Y), {'double', 'cell'})
        error('ERROR: Invalid format for Y.')
    end

%% Define hash functions

    if extflag
        sizeY = Y{2};
    else
        sizeY = size(Y);
        nnzY = nnz(Y);
    end
    N = length(sizeY);

%% Initialize core tensors
    G = cell(N,1); %rng(100)
    for n = 1:N
          G{n} = zeros(R(n) * R(n+1), sizeY(n));
    %     G{n} = rand(R(n), R(n+1)*sizeY(n));    % R=[1 r r ... r 1]; 
    %     G{n} = eye(R(n), R(n+1)*sizeY(n));
    end

%% Main loop: Iterate until convergence, for a maximum of maxiters iterations
    err = []; rel = []; T = [];
    vecY = reshape(double(Y),[],1);
    ns = 1:N;
    YsT = cell(N, 1);
    Gs_hat = cell(N+1,1);  Gs_hat{1} = ones(J1, 1);  Gs_hat{N+1} = ones(1, J1);
    h1int64 = cell(N, 1);    s = cell(N, 1);    
    ite =1;    rel_err = tol +1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for n = 1:N
    %     h1int64{n} = int64(randi(J1, sizeY(n), 1));
    %     s{n} = round(rand(sizeY(n), 1))*2 - 1;
    % end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (ite <= maxiters)  &&  (rel_err > tol)
    
% generate sample matrix
    for n = 1:N
        h1int64{n} = int64(randi(J1, sizeY(n), 1));
        s{n} = round(rand(sizeY(n), 1))*2 - 1;
    end
    
 % generate Yst
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
 % Handle sparse tensor
    if sflag 
        for n = 1:N
            if J1*sizeY(n) < 3*nnzY
                YsT{n} = SparseTensorSketchMatC_git(Y.vals, int64(Y.subs.'), h1int64(ns~=n), s(ns~=n), int64(J1), int64(sizeY(n)), int64(n));
            else
                % This case is for when the sketch dimension is so large that
                % storing YsT{n} in sparse format requires less memory 
                [subs, vals] = Sparse2SparseTensorSketchMatC_git(Y.vals, int64(Y.subs.'), h1int64(ns~=n), s(ns~=n), int64(J1), int64(sizeY(n)), int64(n));
                YsT{n} = sparse(subs(1,:), subs(2,:), vals, J1, sizeY(n));
            end
            %if verbose
                %fprintf('Finished computing sketch %d out of %d...\n', n, N+1)
            %end
        end

 % Use custom external function for computing sketches. Is e.g. used for handling dense large tensor in matfile 
    elseif extflag 
        sketch_func = Y{1};
        if length(Y) > 2
            sketch_params = {J1, h1int64, s, verbose, Y{3:end}};
        else
            sketch_params = {J1, h1int64, s, verbose};
        end
        [YsT] = sketch_func(sketch_params{:});

 % Handle dense normal tensor 
    else 
        for n = 1:N
            YsT{n} = TensorSketchMatC3_git(double(tenmat(Y,n)), h1int64(ns~=n), s(ns~=n), int64(J1)).';
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % right-orthogonalization and generate Gs_hat
    for ind=N:-1:2
        temp1 = reshape(G{ind}, R(ind), R(ind+1)*sizeY(ind));
        [Q, Rfac] = qr(temp1.',0);        
        G{ind-1} = permute(reshape(G{ind-1}, R(ind-1), R(ind), sizeY(ind-1)), [1 3 2]); 
        G{ind-1} = reshape(G{ind-1}, R(ind-1)*sizeY(ind-1), R(ind))*Rfac.';
        R(ind) = size(Q, 2);    
        G{ind-1} = permute(reshape(G{ind-1}, R(ind-1), sizeY(ind-1), R(ind)), [1 3 2]);
        G{ind-1} = reshape(G{ind-1}, R(ind-1)*R(ind), sizeY(ind-1));
        G{ind} = Q.';    G{ind} = reshape(G{ind}, R(ind)*R(ind+1), sizeY(ind));      
        temp2 = fft(countSketch(G{ind}, h1int64{ind}, J1, s{ind}, 1), J1, 2);
        temp2 = reshape(temp2, R(ind), R(ind+1), J1);
        Gs_hat{ind} = zeros(R(ind), J1);
        for j=1:J1
            Gs_hat{ind}(:,j) = temp2(:,:,j)*Gs_hat{ind+1}(:,j);
        end
    end
    
 % inner update
    err_core = 0; 
    tic;
    for ind = 1:N
        H = ifft(dotkron(Gs_hat{ind}, Gs_hat{ind+1}.'));  B = YsT{ind};
        if sigma == 0
           y = pinv((H'*H + sigma*eye(R(ind)*R(ind+1)))) * (H'*B);
%            y = (H'*H) \ (H'*B);
%            y = (H + sigma*eye(size(H))) \ B;
%            y = H \ B;
%            y = (H'*H) \ (H'*B);
%            myfun = @(x) kron(eye(sizeY(ind)), (H'*H)) * x;    vecHtB = H'*B;    y = reshape(pcg(myfun, vecHtB(:)), size(G{ind}));
        else
           y = ((H'*H + sigma*eye(R(ind)*R(ind+1)))) \ (H'*B + sigma*G{ind});
        end
        err_core = max(err_core, norm(y - G{ind},'fro')/norm(y,'fro')); 

 % Update Gs_hat to prepare the next iteration
        if ind < N
            yy = permute(reshape(y, R(ind), R(ind+1), sizeY(ind)), [1 3 2]);
            yy = reshape(yy, R(ind)*sizeY(ind), R(ind+1));
            [Q, Rfac] = qr(yy, 0);     
            G{ind+1} = Rfac*reshape(G{ind+1}, R(ind+1), R(ind+2)*sizeY(ind+1));
            R(ind+1) = size(Q, 2);    G{ind+1} = reshape(G{ind+1}, R(ind+1)*R(ind+2), sizeY(ind+1));
            G{ind} = permute(reshape(Q, R(ind), sizeY(ind), R(ind+1)), [1 3 2]);
            G{ind} = reshape(G{ind}, R(ind)*R(ind+1), sizeY(ind));
            temp = fft(countSketch(G{ind}, h1int64{ind}, J1, s{ind}, 1), J1, 2);
            temp = reshape(temp, R(ind), R(ind+1), J1);
            Gs_hat{ind+1} = zeros(J1, R(ind+1));
            for i=1:J1
                Gs_hat{ind+1}(i,:) = Gs_hat{ind}(i,:)*temp(:,:,i);
            end 
        else
            G{N} = y;
            t = toc; T = [t T];
            vecG = tt_contraction(G);
            rel_err = norm(vecG-vecY)/norm(vecY);
        end
    end
%%%%%%%%%%%%%%%
    fprintf(' Ite %2d: err_core = %7.1e , rel_err = %7.1e \n', ite, err_core, rel_err);
    err = [err err_core];
    rel = [rel rel_err];
    %log_rel = log10(rel);
    ite = ite + 1;
    
end

end
