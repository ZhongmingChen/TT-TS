function vecG = tt_contraction(G)
N=length(G); vecG = 1;

for i=1:N
    R = size(G{i}, 1)/size(vecG, 2); %10/1
    vecG = vecG*reshape(G{i}, size(vecG, 2), []); %1x1000
    vecG = permute(reshape(vecG, size(vecG, 1), R, []), [1 3 2]); %1x10x100
    vecG = reshape(vecG, [], R); %10x100
end

end

% function vec_y = core_contraction(y)
% size_y = size(y);
% N = length(size_y);
% r = zeros(1,N+1); n = zeros(1,N);
% vec_y = y{1};
% for i = 1 : N-1
%     [r(i),n(i),r(i+1)] = size(y{i});
%     vec_y = reshape(vec_y,[],r(i+1))*reshape(y{i+1},r(i+1),[]);
% end
%     vec_y = reshape(vec_y,1,[])';
% end

% error = 0;
% for i = 1000:1000:100000
%     [G, err] = tt_ts(Y, R, i);
%     vecG = tt_contraction(G);
%     rel=norm(vecG'-y,'fro')/norm(y,'fro');
%     error = [error rel];
% end
% figure;plot(x,vecG,'.')