% 3, 11, 17

filenumber = 2;
links = load(sprintf('links/%d.txt',filenumber));
tit_file = fopen(sprintf('titles/%d.txt',filenumber));
titles = textscan(tit_file,'%s');
titles = titles{1};

col_1 = links(:,1);
col_2 = links(:,2);

A = sparse(col_1, col_2, 1);


%% %Task 1

in_degree = ones(1,max(col_1));
out_degree = ones(1,max(col_1));
for i=1:max(col_1)
    out_degree(i) = sum( A(i,:));
    in_degree(i) = sum( A(:,i));
end



[~,B]=maxk(out_degree,5);
[~,C]=maxk(in_degree,5);

fprintf('Out degrees\n')
for i = B
    fprintf('%s ', titles{i})
    fprintf('%.15g , ', out_degree(i)/sum(out_degree))
    fprintf('%.15g\n', in_degree(i)/sum(in_degree))
end
fprintf('\n')
fprintf('In degrees\n')
for i = C
    fprintf('%s ', titles{i})
    fprintf('%.15g, ', in_degree(i)/sum(in_degree))
    fprintf('%.15g\n', out_degree(i)/sum(out_degree))
end

fprintf('\n')

fprintf('%s ', titles{2})
fprintf('%.15g, ', (in_degree(2)+out_degree(2))/(sum(in_degree)+sum(out_degree)))

%B = sort(in_degree, 'descend');
%C = sort(out_degree, 'descend');


%%
%Task 2


mat_hub = transpose(A)*A;
[V, D] = eigs(mat_hub);
[~,idx] = max(diag(D));

hub = abs(V(:,idx));
hub_norm = hub./sum(hub);

[~,F]=maxk(hub_norm,5);

%fprintf('Hub centrality:\n')
for i = F
    fprintf('%s\n', titles{i})
    fprintf('%.15g\n', hub_norm(i))
end

%%

mat_auth = A*transpose(A);
[V, D] = eigs(mat_auth);
[~,idx] = max(diag(D));

auth = abs(V(:,idx));
auth_norm = auth./sum(auth);

[~,L]=maxk(auth_norm,5);

%fprintf('Hub centrality:\n')
for i = L
    fprintf('%s\n', titles{i})
    fprintf('%.15g\n', auth_norm(i))
end

%%
%Task 3

[V, D] = eigs(A);
[lambda_max,idx] = max(diag(D));

alpha = 0.85/abs(lambda_max);

N = max(col_1);
u = ones(N,1);


katz = (1/N)*inv( (eye(N)-alpha*transpose(A)) )*u;

katz_norm = katz./sum(katz);
[~,T]=maxk(katz,5);

for j = T
    fprintf('%s\n', titles{j})
    fprintf('%.15g\n', katz_norm(j))
end

%%
%Task 4

alpha=0.85;
mat_u = u*transpose(u);

H = transpose(A)./out_degree;
P_Rank = alpha*H + ((1-alpha)/N)*mat_u;
P_Rank(isinf(P_Rank)|isnan(P_Rank)) = 1/N; % Replace NaNs and infinite values with zeros

[V, D] = eig(P_Rank);
[lambda_max,idx] = max(diag(D));

vec = abs(V(:,idx));
vec_norm = vec./sum(vec);

[~,F2]=maxk(vec_norm,5);

% PageRank_norm = PageRank./sum(PageRank);
% [~,T]=maxk(PageRank,5);

F2
for j = F2
    fprintf('%s\n', titles{j})
    fprintf('%.15g\n', vec_norm(j))
end


%%
%Task 5



