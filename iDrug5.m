function [S, V, U] = iDrug5(X, w, DV, DS, S, rank1, rank2, para)
% Summary of iDrug function
%   The function solves the cross-domain recommendation for two tasks.
%   Task 1: drug-disease association and Task 2: drug-target assoication.

% INPUT:
% X: contains X{1} for known drug-disease associations matrix, 
%    and X{2} for known drug-target association matrix
% w: the weight constant for defined the weight matrix W
% DV: DV{1} and DV{2} contain the drug-drug similarity for two tasks
% DS: DS{1} contain disease-disease similarity for task 1
%     DS{2} contain target-target similarity for task 2
% S: the mapping matrix for overlap of drugs between two tasks
% rank1, rank2: the rank of the latent factor matrix U{1} and U{2}
% para: contain the regulairzed parameters: [alpha, beta, gamma] 

% OUTPUT:
% U, V: the latent matrix for two tasks:U{1}, U{2} and V{1} V{2}


%%%   Constrcution of the affinity matrix S
for i =1:2
    [DS{i},DV{i},WS{i},WV{i}] = makemanifoldd(X{i});
end

if nargin < 8
    para = [1.00, 1.00, 1.00];
end

for i =1:2
    X{i} = X{i}'
end


alpha = para(1);
beta = para(2);
lambda = para(3);

itermax = 20;

myeps = 1.000000000000000e-09;

l = 500;




%%

S = {};
V = {};

for j = 1:2

    [m,~]=size(X{j});
    S{j}=rand(m,l);
    V{j}=rand(l,m);
    E=X{j}'-(X{j}')*S{j}*V{j};
    a=4*max(diag(E*(E')),myeps);
    U = diag(sqrt(1./a));
    onesm=ones(m,m);
    for i=1:itermax
        %% Update S
        A=X{j}*U*(X{j}');
        S1=A*(V{j}')+(alpha*X{j}*WS{j}*X{j}'+(beta+lambda)*eye(m))*S{j};
        S2=A*S{j}*V{j}*(V{j}')+(alpha*X{j}*DS{j}*(X{j}')+beta*onesm+lambda*S{j}*(S{j}'))*S{j};
        re1=rdivide(S1,S2);
        S{j}=S{j}.*re1;
        %% Update V
        B=(S{j}')*X{j}*U*(X{j}');
        V1=B+alpha*V{j}*WV{j};
        V2=B*S{j}*V{j}+alpha*V{j}*DV{j};
        re2=rdivide(V1,V2);
        V{j}=V{j}.*re2;
        %% Update U
        E=X{j}'-(X{j}')*S{j}*V{j};
        a=4*max(diag(E*(E')),myeps);
        U = diag(sqrt(1./a));
    end
    tempVector = sum(S{j}.^2, 2);

end

%%
% 
% % Initilization
% rank = [rank1, rank2];
% 
% for i = 1:2
%     [row, col] = size(X{i});
%     U{i} = rand(row, rank(i));
%     U{i} = U{i}/sqrt(trace(U{i}'*U{i}));
%     V{i} = rand(col, rank(i));
%     V{i} = V{i}/sqrt(trace(V{i}'*V{i}));
% end
% 
% % compute R{i}
% for i = 1:2
%     R{i} = U{i}*V{i}';
%     R{i}(X{i}==0) = 0;
% end
% 
% Iter = 1;
% maxIter = 200;
% epsilon = 1e-4;
% relErr = 10e8;
% Jopt1 = objectiveValue(X, U, V, w, DV, DS, S, para);
% objs = [Jopt1];
% 
% 
% % Update U{i} and V{i} alternately
% SS = S'*S;
% while Iter <= maxIter && relErr > epsilon
%     
%     % Update U{1}
%     [row1, ~] = size(DV{1});
%     num = X{1}*V{1} + para(1)*DV{1}*U{1};
%     num = num + 2*para(2)*S'*U{2}*U{2}'*S*U{1};
%     den = ((1-w^2)*R{1} + w^2*U{1}*V{1}')*V{1} + para(1)*spdiags(sum(DV{1},2),0,row1,row1)*U{1};
%     den = den + 2*para(2)*S'*S*U{1}*U{1}'*SS*U{1} + 0.5*para(3);
%     U{1} = U{1}.*((num./den).^(0.5));
%     
%     % Update U{2}
%     [row2, ~] = size(DV{2});
%     num = X{2}*V{2} + para(1)*DV{2}*U{2};
%     num = num + 2*para(2)*S*U{1}*U{1}'*S'*U{2};
%     den = ((1-w^2)*R{2} + w^2*U{2}*V{2}')*V{2} + para(1)*spdiags(sum(DV{2},2),0,row2,row2)*U{2};
%     den = den + 2*para(2)*U{2}*U{2}'*U{2} + 0.5*para(3);
%     U{2} = U{2}.*((num./den).^(0.5));
%     
%     % Update V{i}
%     for i = 1:2
%          [row, ~] = size(DS{i});
%          num = X{i}'*U{i} + para(1)*DS{i}*V{i};
%          den = ((1-w^2)*R{i}' + w^2*V{i}*U{i}')*U{i} +  para(1)*spdiags(sum(DS{i},2),0,row,row)*V{i} + 0.5*para(3);
%          V{i} = V{i}.*((num./den).^(0.5));
%     end
% 
%     for i = 1:2
%         R{i} = U{i}*V{i}';
%         R{i}(X{i}==0) = 0;
%     end
%     
%     Jopt2 = objectiveValue(X, U, V, w, DV, DS, S, para);
%     relErr = Jopt1 - Jopt2;
%     objs = [objs, Jopt2];
%     Jopt1 = Jopt2;
%     Iter = Iter + 1;
% 
% end
% 
% 
% end
%       
% 
% %%% compute the value of the objective function
% function J = objectiveValue(X, U, V, w, DV, DS, S, para)
% % construct the W weight matrix
% for i = 1:2
%     W{i} = ones(size(X{i})) * w;
%     W{i}(X{i} > 0) = 1;
% end
% 
% % construct the laplacian matrix WV and WS
% for i = 1:2
%     [row1, ~] = size(DV{i});
%     [row2, ~] = size(DS{i});
%     WV{i} = -DV{i} + spdiags(sum(DV{i},2),0,row1,row1);
%     WS{i} = -DS{i} + spdiags(sum(DS{i},2),0,row2,row2);
% end
% 
% 
% J = 0;
% 
% for i = 1:2
%     J = J + norm(W{i}.*(X{i} - U{i}*V{i}'), 'fro')^2;
%     J = J + para(1)* (trace(U{i}'* WV{i} * U{i}) + trace(V{i}'* WS{i} * V{i}));
%     J = J + para(3) * (sum(U{i}(:)) + sum(V{i}(:)));
% end
% 
% J = J + para(2) * norm(S*U{1}*(S*U{1})' - U{2}*U{2}', 'fro')^2;
% 
% end
