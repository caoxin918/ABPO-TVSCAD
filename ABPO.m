function[x_k] = PSB(Nodes,b_analytic,G)
NumNode = size(G,2);
G_inv = G';
G_inv_G = G_inv * G;
% ABPO
E_L2_sum = zeros(NumNode,1);
E_cos_sum = zeros(NumNode,1);
E_L2_0 = 1000;
E_cos_0 = 0;
% 初始化x
x_k = pinv(G)*b_analytic;
y_k = zeros(NumNode-1,1);
z_k = zeros(NumNode-1,1);
w = zeros(NumNode-1,1);

D = zeros(NumNode - 1, NumNode);


% 给每个中心顶点找k个最近的结点
% for i = 1:NumNode-1
%     [Idx,~] = knnsearch(Nodes, Nodes(i,:),'K',4);
%     D(i,i) = -1;
%     D(i,Idx) = 1/3;
% end

% 一阶差分矩阵
D = spdiags( [ -ones( NumNode, 1 ) ones( NumNode, 1 ) ], 0 : 1, NumNode - 1 , NumNode );
D_inv = D';
D_inv_D = D_inv * D;
% 迭代次数
iter = 1;
% 算法中的一些参数
lambda = 0.05;
eta = 0.5;
alpha = 3.7;
% 
while 1
    gamma1 = 0.03;
    gamma2 = alpha * gamma1;
%     gamma2 = max(0.2*0.85^(k-1),0.1);
    
%     z = Dx,是x向量的差分
    z = D*x_k;
    TEMP = G_inv * b_analytic + eta * D_inv * z - D_inv * y_k;
    %     开始算法
    x_0 = x_k;
    x_k = (G_inv_G + eta * D_inv_D) \ TEMP;
    D_x_k = D * x_k;
    t = D_x_k + y_k / eta;
    
    left_idx = find(abs(t) <= 2 * gamma1);
    mid_idx = find(abs(t) > 2 * gamma1 & abs(t) <= gamma2);
    right_idx = find(abs(t) > gamma2);
    
    z_k(left_idx) = sign(t(left_idx)) .* max(abs(t(left_idx)) - gamma1,0);
    z_k(mid_idx) = ( (alpha - 1) * t(mid_idx) - sign(t(mid_idx)) * alpha * gamma1 ) / (alpha - 2);
    z_k(right_idx) = t(right_idx);

    y_k = y_k - eta * (D_x_k - z_k);
    
    iter = iter + 1;
    x_k(find(x_k < 0)) = 0;
    
    E_L2 = norm(G * x_k - b_analytic);
    E_cos = sum(G * x_k .* b_analytic) / (norm(G * x_k) * norm(b_analytic));
    
    E_L2_sum = E_L2_sum + E_L2;
    E_cos_sum = E_cos_sum + E_cos;
    
    P_L2 = (1 ./ E_L2_sum) / sum(1 ./ E_L2_sum);
    P_cos = E_cos_sum / sum(E_cos_sum);
    P_err = (P_L2 + P_cos) / 2;
    P_X = x_k .* P_err / sum(x_k .* P_err);
    
    
    if E_cos < E_cos_0 || E_L2 > E_L2_0
        break;
    end
    E_cos_0 = E_cos;
    E_L2_0 = E_L2;
    x = x_k;
end