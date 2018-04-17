% c00194582 caizhongheng 2018.04.07
clear
clc
close all
rng(257);

%% setting aera
bit_grp = 1000;
k = 6;
m = 9;
n = k + m;

Hb = [1 0 1 -1 -1; -1 0 0 0 -1;0 -1 -1 0 0];
Z  = eye(3);
H  = zeros(size(Hb,1)*length(Z),size(Hb,2)*length(Z));

for idx_row = 1:size(Hb,1)
    for idx_col = 1:size(Hb,2)
        if(Hb(idx_row,idx_col)==-1)
            H((idx_row-1)*length(Z)+(1:length(Z)),(idx_col-1)*length(Z)+(1:length(Z))) = zeros(size(Z));
        else
            H((idx_row-1)*length(Z)+(1:length(Z)),(idx_col-1)*length(Z)+(1:length(Z))) = circshift(Z,[0 Hb(idx_row,idx_col)]);
        end
    end
end

HI = H(:,1:k);
HP = H(:,k+1:end);

%% create the bit data in
bit_in = randi([0 1],k*bit_grp,1);
bit_in = reshape(bit_in,k,[]);

%% use LDPC to encode
u = mod(HI*bit_in,2);
p = mod(mod(inv(HP),2)*u,2);
% check code
% H_check = mod(HI*bit_in + HP*p,2);

enc_out = [bit_in;p];

%% mapping and calc llr
llr_data = 2*double(~enc_out)-1;
% 因为llr=ln(APP(X=0)/APP(X=1))，这里为了简化llr的计算，所以对enc_out做了一个取反动作，然后直接映射成BPSK信号，作为LLR值。
%% add noise
noise = randn(size(llr_data))+1j*randn(size(llr_data));
noise = noise/2.4;% 2.9是调整噪声功率的比值，目的是让BPSK信号既出现纠前误码又不至于超出纠前门限
tmp_data = llr_data + noise;

scatterplot(tmp_data(:))
llr_data = real(llr_data + noise);

%% decode

iter_num = 5; % 迭代次数
dec_grp = bit_grp;


% 找出H矩阵中的非0元素，方便迭代用。
[node_row,node_col] = find(H==1);
node_array = [node_row node_col];
clear node_row node_col
% C2V表示校验方程。
C2V_slct     = zeros(size(H,1),length(find(node_array(:,1)==1)));
for idx=1:size(C2V_slct,1)
    % 找到H矩阵每行中为1的位置
    Q_adress      = find(node_array(:,1)==idx);
    % 将该位置存入Q_slct矩阵中
    C2V_slct(idx,:) = node_array(Q_adress,2);
end

% V2C，行为变量节点，列为和该变量节点相连结的校验节点
V2C_slct     = zeros(size(H,2),2);
for idx=1:size(V2C_slct,1)
    Q_adress = find(node_array(:,2)==idx);
    V2C_slct(idx,:) = node_array(Q_adress,1);
end





llr_data_in = llr_data(:,1:dec_grp);% first col of data
dec_out     = zeros(n,dec_grp);% 硬判输出结果
llr_array = zeros(size(llr_data_in,1),iter_num);

% keyboard;
% llr_data_in(end) = 1.2;

for grp_idx=1:dec_grp
    tmp_llr_data = llr_data_in(:,grp_idx);
    
    V = zeros(15,2);
    APP = tmp_llr_data;
    C = zeros(9,3);
    tmp_dec_out = zeros(15,1);
    
    for iter_idx=1:iter_num
        % 更新V的概率信息
        for idx=1:n
            [R_row,R_col] = find(C2V_slct==idx);
            C_slct        = [R_row,R_col];
            %             V(idx) = APP(idx);
            if(idx<=12)
                V(idx,1) = APP(idx) - C(C_slct(1,1),C_slct(1,2));
                V(idx,2) = APP(idx) - C(C_slct(2,1),C_slct(2,2));
            else
                V(idx,1) = APP(idx) - C(C_slct(1,1),C_slct(1,2));
            end
        end
        
        % 校验节点信息更新：Z1到Z9
        for idx=1:m
            %             if(idx==m)
            %                 keyboard;
            %             end
            idx_array_1 = V2C_slct(C2V_slct(idx,1),:);
            new_idx_1 = find(idx_array_1==idx,1);
            idx_array_2 = V2C_slct(C2V_slct(idx,2),:);
            new_idx_2 = find(idx_array_2==idx,1);
            idx_array_3 = V2C_slct(C2V_slct(idx,3),:);
            new_idx_3 = find(idx_array_3==idx,1);
            
            C(idx,1) = sign(V(C2V_slct(idx,2),new_idx_2))*sign(V(C2V_slct(idx,3),new_idx_3))*min(abs(V(C2V_slct(idx,2),new_idx_2)),abs(V(C2V_slct(idx,3),new_idx_3)));
            C(idx,2) = sign(V(C2V_slct(idx,1),new_idx_1))*sign(V(C2V_slct(idx,3),new_idx_3))*min(abs(V(C2V_slct(idx,1),new_idx_1)),abs(V(C2V_slct(idx,3),new_idx_3)));
            C(idx,3) = sign(V(C2V_slct(idx,1),new_idx_1))*sign(V(C2V_slct(idx,2),new_idx_2))*min(abs(V(C2V_slct(idx,1),new_idx_1)),abs(V(C2V_slct(idx,2),new_idx_2)));

            %             C(idx,1) = 2*atanh(tanh(V(C2V_slct(idx,2))/2)*tanh(V(C2V_slct(idx,3))/2));
            %             C(idx,2) = 2*atanh(tanh(V(C2V_slct(idx,1))/2)*tanh(V(C2V_slct(idx,3))/2));
            %             C(idx,3) = 2*atanh(tanh(V(C2V_slct(idx,1))/2)*tanh(V(C2V_slct(idx,2))/2));
            %             R_1(idx,1) = sign(V(C2V_slct(idx,2)))*sign(V(C2V_slct(idx,3)))*min(abs(V(C2V_slct(idx,2))),abs(V(C2V_slct(idx,3))));
            %             R_1(idx,2) = sign(V(C2V_slct(idx,1)))*sign(V(C2V_slct(idx,3)))*min(abs(V(C2V_slct(idx,1))),abs(V(C2V_slct(idx,3))));
            %             R_1(idx,3) = sign(V(C2V_slct(idx,1)))*sign(V(C2V_slct(idx,2)))*min(abs(V(C2V_slct(idx,2))),abs(V(C2V_slct(idx,2))));
        end
        
        % 变量节点信息更新：P1到P15
        for idx=1:n
            [R_row,R_col] = find(C2V_slct==idx);
            C_slct        = [R_row,R_col];
            if(idx<=12)
                APP(idx) = APP(idx) + C(C_slct(1,1),C_slct(1,2)) + C(C_slct(2,1),C_slct(2,2));
            else
                APP(idx) = APP(idx) + C(C_slct(1,1),C_slct(1,2));
            end
            % 硬判
            if(APP(idx)>=0)
                tmp_dec_out(idx) = 0;
            else
                tmp_dec_out(idx) = 1;
            end
        end
        llr_array(:,iter_idx) = APP;
        check_out = mod(H*tmp_dec_out,2);
        check_out_sum = sum(check_out);
        %         if(check_out_sum~=0)
        %             keyboard
        %         else
        %         end
        
    end
    dec_out(:,grp_idx) = tmp_dec_out;
end

check = abs(dec_out(1:15,:) - enc_out(1:15,1:dec_grp));
sum(sum(check))
