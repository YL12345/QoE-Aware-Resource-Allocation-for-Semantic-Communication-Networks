function [QoE_sum,state_u]=SCMatching(BS_position,SINR_single,DeepSC_table,SINR_Bi,VQA_table,N,radius,N_cell,N_channels,shadow_factor,Nr,P_range,P_noise,I_th,H_S,K_S,bandwidth,H_Bi_text,H_Bi_image,K_Bi_text,K_Bi_image,G_th)
% solve the sum QoE of the proposed matching method; single cell optimization
% Input: BS_position: cell:N_cell*1, the positions of BSs
%        SINR_single: the sinr range of single-modal users
%        DeepSC_table: the performance table of single-modal users; row: number of symbols; column:snr
%        SINR_Bi: all possible combinations of SINR of two users
%        VQA_table: VQA_table=cell(length(K_Bi_image),length(K_Bi_text)); in each cell, row:snr of text user; column:snr of image user
%        N: N_cell*3, number of users; the first column: N_S; the second column: N_Bi; the third column: N_D=N_S+N_Bi
%        radius: the radius of each cell
%        N_cell: number of cells
%        N_channels: number of channels
%        shadow_factor: shadow factor for large scale fading
%        Nr: number of receive antennas
%        P_range: the range of transmit power (mW)
%        P_noise: the noise power (mW)
%        I_th: the interference threshold from the current cell to other
%        H_S: semantic entropy of single-modal user
%        K_S: possible values of of number of semantic symbols for DeepSC
%        bandwidth: bandwidth of each channel
%        H_Bi_text/H_Bi_image: semantic entropy of bimodal user
%        K_Bi_text/K_Bi_image: possible values of of number of semantic symbols for VQA
%       G_th=0.5; %the minimum score of phi and si for all users
% Output: QoE_sum: sum QoE of all users in all cells
% tic
QoE_sum=0;
P_sum=0;

w_phi=cell(N_cell,1); % random weight of phi for every user; w_si=1-w_phi;
para_S=cell(N_cell,1);
para_Bi_text=cell(N_cell,1);
para_Bi_image=cell(N_cell,1);


%% locations of users
D_position=cell(N_cell,1); %the positions of users in N_cell cells
for n_cell=1:1:N_cell
    N_D=N(n_cell,3); % number of all users in cell i
    D_position{n_cell}=zeros(N_D,2); %存储BS1下的用户位置坐标
%     radius_d=(radius-50)+50*sqrt(rand(N_D,1)); %radius of users position
    radius_d=radius*sqrt(rand(N_D,1)); %radius of users position
    phase=rand(N_D,1)*2*pi; %phases of users
    D_position{n_cell}(:,1)=radius_d.*cos(phase)+BS_position{n_cell}(1); %x coodinate of users
    D_position{n_cell}(:,2)=radius_d.*sin(phase)+BS_position{n_cell}(2); %y coodinate of users
end

%% channel: large scale fading, including pathloss and shadowing; small scale fading: Rayleigh fading;
h=cell(N_cell,N_channels); %channel fading coefficient, not channel gains, including large scale fading and small scale fading
h_large=cell(N_cell,1); %large scale fading of users in different cells
h_small=cell(N_cell,N_channels); %small scale fading of users in different cells over diffrent channels
for n_cell=1:1:N_cell %the cell where users are deployed
    N_D=N(n_cell,3); % number of all users in cell i
    h_large{n_cell}=zeros(N_D,N_cell); %large scale fading from users in n_cell to every BSs
    for n_cell_1=1:1:N_cell %the cell of BSs
        d=sqrt((D_position{n_cell}(:,1)-BS_position{n_cell_1}(1)).^2 + (D_position{n_cell}(:,2)-BS_position{n_cell_1}(2)).^2);%ditance from users in n_cell to BS in n_cell_1
        pl=128.1 + 37.6 * log10(d/1000); %pathloss
        h_large{n_cell}(:,n_cell_1)=10.^(-(pl + shadow_factor)/10);
    end
    for n_channels=1:1:N_channels %for each channel
       h_small{n_cell,n_channels}=zeros(N_cell,Nr,N_D); %from users to different BSs; for each receive antenna
       h{n_cell,n_channels}=zeros(N_cell,Nr,N_D);
       for n_d=1:1:N_D %for each user in n_cell, get the small scale fading coefficients
           h_real = randn(N_cell,Nr)/sqrt(2);
           h_image = randn(N_cell,Nr)/sqrt(2);
           h_small{n_cell,n_channels}(:,:,n_d)=h_real+h_image*1i;
           for n_cell_1=1:1:N_cell
               h{n_cell,n_channels}(n_cell_1,:,n_d)=h_small{n_cell,n_channels}(n_cell_1,:,n_d)*sqrt(h_large{n_cell}(n_d,n_cell_1)); %multipy small scale fading by corresponding large scale fading
           end
       end
    end
end

% obtain h for single cell optimization
h_SC=cell(N_cell,N_channels); %channel fading coefficient, not channel gains, including large scale fading and small scale fading
for n_cell=1:1:N_cell
    for n_channels=1:1:N_channels
        h_SC{n_cell,n_channels}=h{n_cell,n_channels}(n_cell,:,:);
    end
end

%% the maximum QoE and optimal k under every possible SNR for each user group
k=cell(N_cell,1); % the optimal k of each single-modal user under SINR_single
QoE_k=cell(N_cell,1); % the optimal k of each single-modal user under SINR_single

%% initialize the states of user groups and channels
state_u=cell(N_cell,1); % states of users of N_cell cells
state_c=cell(N_cell,1); % states of channels of N_cell cells
for n_cell=1:1:N_cell
    N_D=N(n_cell,3); % number of all users in cell i
    N_S=N(n_cell,1)+max(N_D,N_channels)-N_D; % number of single-modal users; add virtual users
    N_Bi=N(n_cell,2); % number of bimodal users in cell i
    state_u{n_cell}=cell(N_Bi/2+N_S,1); % include N_Bi/2+N_S user groups
    state_c{n_cell}=zeros(max(N_D,N_channels),2); % add virtual channels; row: channels; 1st column: the current QoE; 2rd column: index of user group
end

%% solve the optimal resource allocation for each cell
for n_cell=1:1:N_cell
    [QoE,state_u(n_cell),state_c(n_cell),QoE_k(n_cell),k(n_cell),w_phi(n_cell),para_S(n_cell),para_Bi_text(n_cell),para_Bi_image(n_cell)]=SCMinPowerMatching(h_SC(n_cell,:),SINR_single,DeepSC_table,SINR_Bi,VQA_table,N(n_cell,:),1,N_channels,P_range,P_noise,H_S,K_S,bandwidth,H_Bi_text,H_Bi_image,K_Bi_text,K_Bi_image,G_th);
end

%% solve the actual QoE
% save the optimization results in this single cell optimization method,so as to calculate the final QoE considering the inter-cell interference
power_result=cell(N_cell,1);
channel_result=cell(N_cell,1);
k_result=cell(N_cell,1);
for n_cell=1:1:N_cell
    N_S=N(n_cell,1);
    N_Bi=N(n_cell,2);
    power_result{n_cell}=zeros(1,N_S+N_Bi);
    channel_result{n_cell}=zeros(1,N_S+N_Bi);
    k_result{n_cell}=zeros(1,N_S+N_Bi);
    for q=1:1:N_Bi/2 % for bimodal user
        if state_u{n_cell}{q}(4,1)==0 % QoE = 0; no channel assigned
           k_result{n_cell}(2*q-1:2*q)=[0,0];
           channel_result{n_cell}(2*q-1:2*q)=[0,0];
           power_result{n_cell}(2*q-1:2*q)=[0,0];
           continue
        end
       sinr=state_u{n_cell}{q}(3,:); 
       sinr_temp=sinr.'-SINR_Bi;
       temp=max(find(sinr_temp(1,:)>=0&sinr_temp(2,:)>=0)); % the index of SINR in SINR_Bi
       k_result{n_cell}(2*q-1:2*q)=k{n_cell}{q}(:,temp)';
       channel_result{n_cell}(2*q-1:2*q)=state_u{n_cell}{q}(1,:);
       power_result{n_cell}(2*q-1:2*q)=state_u{n_cell}{q}(2,:);
    end
    for n_s=N_Bi/2+1:1:N_Bi/2+N_S % for each single-modal user
       if state_u{n_cell}{n_s}(4,1)==0 % QoE = 0; no channel assigned
           k_result{n_cell}(n_s+N_Bi/2)=0;
           channel_result{n_cell}(n_s+N_Bi/2)=0;
           power_result{n_cell}(n_s+N_Bi/2)=0;
           continue
        end
       sinr=state_u{n_cell}{n_s}(3);
       sinr_temp=sinr-SINR_single;
       temp=max(find(sinr_temp>=0)); % the index of SINR in SINR_Bi
       k_result{n_cell}(n_s+N_Bi/2)=k{n_cell}{n_s}(temp);
       channel_result{n_cell}(n_s+N_Bi/2)=state_u{n_cell}{n_s}(1);
       power_result{n_cell}(n_s+N_Bi/2)=state_u{n_cell}{n_s}(2);
    end
end

% calculate the actual QoE
[QoE_real]=QoEReal(K_Bi_text, K_Bi_image, para_Bi_text, para_Bi_image, H_Bi_text, H_Bi_image, G_th, w_phi, para_S, bandwidth, H_S, N, channel_result, k_result, power_result, N_cell, P_noise, h);


%% save the results
for n_cell=1:1:N_cell
    QoE_sum=QoE_sum+sum(QoE_real{n_cell});
end

% toc
end