% 3, 11, 17

filenumber = 2;
links = load(sprintf('links/%d.txt',filenumber));
tit_file = fopen(sprintf('titles/%d.txt',filenumber));
titles = textscan(tit_file,'%s');
titles = titles{1};

col_1 = links(:,1);
col_2 = links(:,2);

A = transpose( sparse(col_1, col_2, 1, 3000, 3000) );


%% %Task 1

in_degree = ones(1,max(col_2));
out_degree = ones(1,max(col_2));
for i=1:max(col_1)
    in_degree(i) = sum( A(i,:));
    out_degree(i) = sum( A(:,i));
end

[~,C]=maxk(out_degree,5);
[~,B]=maxk(in_degree,5);

fprintf('Top In degrees:\n')
for i = B
    fprintf('%s ', titles{i})
    fprintf('%.6g  ,  ', in_degree(i)/sum(in_degree))
    fprintf('%.6g\n', out_degree(i)/sum(out_degree))
    
end
fprintf('\n')
fprintf('Top out degrees\n')
for i = C
    fprintf('%s ', titles{i})
    fprintf('%.6g  ,  ', in_degree(i)/sum(in_degree))
    fprintf('%.6g\n', out_degree(i)/sum(out_degree))
    
end
fprintf('\n')


%%
%Task 2

[V1, D1] = eigs( transpose(A)*A );
[~,id1] = max(diag(D1));

hub = abs(V1(:,id1));
hub_norm = hub./sum(hub);

[~,F]=maxk(hub_norm,5);


[V2, D2] = eigs( A*transpose(A) );
[~,id2] = max(diag(D2));

auth = abs(V2(:,id2));
auth_norm = auth./sum(auth);

[~,L]=maxk(auth_norm,5);

fprintf('Top Hubs:\n')
for i = F
    fprintf('%s\n', titles{i})
    fprintf('%.5g\n', hub_norm(i))
    fprintf('%.5g\n', auth_norm(i))
end
fprintf('\n')
fprintf('Top authorities:\n')
for i = L
    fprintf('%s\n', titles{i})
    fprintf('%.5g\n', auth_norm(i))
    fprintf('%.5g\n', hub_norm(i))
end



%%
%Task 3

[V3, D3] = eigs(A);
[lambda_max,id3] = max(diag(D3));

alpha = 0.85/abs(lambda_max);

N = max(col_2);
u = ones(N,1);

katz = (1/N)*inv( (eye(N)-alpha*A) )*u;

katz_norm = katz./sum(katz);
[~,T]=maxk(katz,5);

fprintf('Top Katz:\n')
for j = T
    fprintf('%s\n', titles{j})
    fprintf('%.5g\n', katz_norm(j))
end

%%
%Task 4

alpha=0.85;
mat_u = u*transpose(u);

H = A./out_degree;
P_Rank = alpha*H + ((1-alpha)/N)*mat_u;
P_Rank(isinf(P_Rank)|isnan(P_Rank)) = 1/N; % Replace NaNs and infinite values with zeros

[V4, D4] = eig(P_Rank);
[lambda_max,id4] = max(diag(D4));

vec = abs(V4(:,id4));
vec_norm = vec./sum(vec);

[~,F]=maxk(vec_norm,5);

fprintf(' aplha = 0.85, Top Page Rank:\n')
for j = F
    fprintf('%s\n', titles{j})
    fprintf('%.4g\n', vec_norm(j))
end

fprintf('\n')

N = max(col_2);
H(isinf(H)|isnan(H)) = 1/N; % Replace NaNs and infinite values with zeros

page_rank = (1/N)*inv( (eye(N) - alpha*H) )*u;

page_rank_norm = page_rank./sum(page_rank);

[~,K]=maxk(page_rank_norm,5);

top_three = [K(1), K(2), K(3)];
for j = K
    fprintf('%s\n', titles{j})
    fprintf('%.4g\n', page_rank_norm(j))
end

%%
%Task 5

alpha=0.85;
N = max(col_2);

iterations = 100;

top1 = ones(1,iterations);
top2 = ones(1,iterations);
top3 = ones(1,iterations);

pRank = ones(N,1);

mat_u = u*transpose(u);
G = alpha*H + ( (1-alpha)/N )*mat_u;

for k = 1:iterations
    
    pRank = G*pRank;
    
    pRank_norm = pRank./sum(pRank);

    top1(k) = pRank_norm(top_three(1));
    top2(k) = pRank_norm(top_three(2));
    top3(k) = pRank_norm(top_three(3));

    [~,W]=maxk(pRank_norm,5);

    fprintf('Top pageRank itteration %.1g:\n', k)
    for j = W
        fprintf('%s', titles{j})
        fprintf(' :  %.4g\n', pRank_norm(j))
    end
    fprintf('\n')
end
%Task 6 plot

figure(1)
plot(1:k, top1)
hold on
plot( zeros(k,1) + page_rank_norm( top_three(1) ) )
title(titles{top_three(1)})

figure(2)
plot(1:k, top2)
hold on
plot( zeros(k,1) + page_rank_norm( top_three(2) ) )
title(titles{top_three(2)})

figure(3)
plot(1:k, top3)
hold on
plot( zeros(k,1) + page_rank_norm( top_three(3) ) )
title(titles{top_three(3)})
