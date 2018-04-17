% c00194582 caizhongheng 2018.04.14
% add layer_decoding and modify the encoding process

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
HIb = Hb(:,1:2);

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

u = zeros(size(Hb,1)*length(Z),bit_grp);
%% use LDPC to encode
for grp_idx = 1:size(bit_in,2)
    Mb = reshape(bit_in(:,grp_idx),length(Z),size(bit_in,1)/length(Z));
    tmp_u = zeros(size(Hb,1)*length(Z),2);
    for row_idx = 1:size(Mb,1)
        for col_idx=1:size(Mb,2)
            if(HIb(row_idx,col_idx)==-1)
                tmp_u((1:length(Z))+length(Z)*(row_idx-1),col_idx) = zeros(length(Z),1);
            else
                tmp_u((1:length(Z))+length(Z)*(row_idx-1),col_idx) = circshift(Mb(:,col_idx),-HIb(row_idx,col_idx));
            end
        end
    end
    tmp_u = mod(sum(tmp_u,2),2);
    u(:,grp_idx) = tmp_u;
end

% u = mod(HI*bit_in,2); u的计算过程也可以简写为这个式子


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
dec_grp  = bit_grp;
llr_data_in = llr_data(:,1:dec_grp);% first col of data
dec_out     = zeros(n,dec_grp);% 硬判输出结果
llr_array   = zeros(size(llr_data_in,1),iter_num);

% keyboard;
% llr_data_in(end) = 1.2;

for grp_idx=1:dec_grp
    tmp_llr_data = llr_data_in(:,grp_idx);
    V = zeros(size(Hb,1),3); %列为校验矩阵的最大行重3
    APP = reshape(tmp_llr_data,length(Z),size(tmp_llr_data,1)/length(Z));
    tmp_dec_out = zeros(15,1);
    
    for iter_idx=1:iter_num
        C = zeros(size(Hb,1),3);
        C_total = repmat(C,length(Z),1);%每次迭代独立进行，所以需要先将C_total清零
        % 更新V的概率信息
        for row_idx = 1:size(Hb,1)% 按Hb逐行更新
            V = zeros(size(Hb,1),3); %列为校验矩阵的最大行重3，因为是不规则的ldpc校验矩阵，所以在每次使用V之前都要先清0一把
            V_idx = 0;
            %% 按Hb的行将APP的值进行循环移位，送给V矩阵
            for col_idx=1:size(Hb,2)
                if(Hb(row_idx,col_idx)~=-1)
                    V_idx = V_idx + 1;
                    V(:,V_idx) = circshift(APP(:,col_idx),-Hb(row_idx,col_idx));
                    % 在V矩阵里面减去上一次迭代送过来的C矩阵
                    V(:,V_idx) = V(:,V_idx) - C_total((1:length(Z))+length(Z)*(row_idx-1),V_idx);
                else
                end
            end
            %% 按Hb的行计算C矩阵
            C(:,1) = sign(V(:,2)).*sign(V(:,3)).*min(abs(V(:,2)),abs(V(:,3)));
            C(:,2) = sign(V(:,1)).*sign(V(:,3)).*min(abs(V(:,1)),abs(V(:,3)));
            C(:,3) = sign(V(:,1)).*sign(V(:,2)).*min(abs(V(:,1)),abs(V(:,2)));
            
            C_total((1:length(Z))+length(Z)*(row_idx-1),:) = C;
            
            %% 将计算出来的C矩阵回加到V矩阵上
            tmp_V = V + C;
            %% 在V矩阵上按照Hb的行进行循环反移位
            % 找出该row_idx下Hb不为-1的列数。
            return_col_idx = find(Hb(row_idx,:)~=-1);
            V_idx = 0;
            for col_idx=return_col_idx
                V_idx = V_idx + 1;
                V(:,V_idx) = circshift(tmp_V(:,V_idx),Hb(row_idx,col_idx));
            end
            
            %% 将更新后的V放回llr的ram APP上，对应回填
            APP(:,return_col_idx) = V;
        end
        
        serial_APP = APP(:);
        for idx=1:length(serial_APP)
            if(serial_APP(idx)>=0)
                tmp_dec_out(idx) = 0;
            else
                tmp_dec_out(idx) = 1;
            end
        end
        llr_array(:,iter_idx) = serial_APP;
    end
    
    %% 按照Hb逐行更新完成后，对APP进行硬判
    
    check_out = mod(H*tmp_dec_out,2);
    check_out_sum = sum(check_out);
    dec_out(:,grp_idx) = tmp_dec_out;
    
end
check = abs(dec_out(1:15,:) - enc_out(1:15,1:dec_grp));
sum(sum(check))





