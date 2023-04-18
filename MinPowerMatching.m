function [QoE_sum,state_u,P_sum]=MinPowerMatching(BS_position,SINR_single,DeepSC_table,SINR_Bi,VQA_table,N,radius,N_cell,N_channels,shadow_factor,Nr,P_range,P_noise,I_th,H_S,K_S,bandwidth,H_Bi_text,H_Bi_image,K_Bi_text,K_Bi_image,G_th)
% solve the sum QoE of the proposed matching method
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
%         state_u: the states of all user groups
%         P_sum: sum power assumption of all users
% tic
QoE_sum=0;
P_sum=0;

%% QoE related parameters for each user
w_phi=cell(N_cell,1); % save the random weights of phi for every user; w_si=1-w_phi;
para_S=cell(N_cell,1); % save parameters for the single-modal user
para_Bi_text=cell(N_cell,1); % parameters for the text transmission user in bimodal task
para_Bi_image=cell(N_cell,1); % parameters for the image transmission user in bimodal task
for i=1:1:N_cell
    N_D=N(i,3); % number of all users in cell i
    N_S=N(i,1); % number of single-modal users in cell i
    N_Bi=N(i,2); % number of bimodal users in cell i
    w_phi{i}=rand(N_D,1);% random weight of phi for every user; w_si=1-w_phi;
    %text transmisstion: 
    %semantic emit rate is a random value from 50 to 70 in Ksuts/s;
    %beta follows CN(0.2,0.05^2)
    phi_S=50+(70-50)*rand(N_S,1);
    beta_S=normrnd(0.2,0.05,N_S,1);
    phi_Btext=50+(70-50)*rand(N_Bi/2,1);
    beta_Btext=normrnd(0.2,0.05,N_Bi/2,1);
    %image transmission:
    %semantic emit rate is a random value from 80 to 100 in Ksuts/s; 
    %beta follows CN(0.1,0.02^2)
    phi_Bimage=80+(100-80)*rand(N_Bi/2,1);
    beta_Bimage=normrnd(0.1,0.02,N_Bi/2,1);
    %si is a random value from 0.8 to 0.9;
    %lamda follows CN(55,2.5^2)
    si=0.8+(0.9-0.8)*rand(N_D,1); %the si threshold of users in a bimodal pair should vary with each other 
    lamda=normrnd(55,2.5,N_D,1); %

    para_S{i}=[beta_S,phi_S,lamda(1:N_S,:),si(1:N_S,:)]; %columns:beta_S;phi_S;lamda_S;si_S
    para_Bi_text{i}=[beta_Btext,phi_Btext,lamda((N_S+1):(N_S+N_Bi/2),:),si((N_S+1):(N_S+N_Bi/2),:)];%colums:beta_Btext;phi_Btext;lamda_Btext;si_Btext
    para_Bi_image{i}=[beta_Bimage,phi_Bimage,lamda((N_S+N_Bi/2+1):N_D,:),si((N_S+N_Bi/2+1):N_D,:)];%colums:beta_Bimage;phi_Bimage;lamda_Bimage;si_Bimage
end

%% locations of users
D_position=cell(N_cell,1); %the positions of users in N_cell cells
for n_cell=1:1:N_cell
    N_D=N(n_cell,3); % number of all users in cell i
    D_position{n_cell}=zeros(N_D,2); % the position of users in n_cell
    radius_d=radius*sqrt(rand(N_D,1)); % radius of users position
    phase=rand(N_D,1)*2*pi; % phases of users
    D_position{n_cell}(:,1)=radius_d.*cos(phase)+BS_position{n_cell}(1); % x coodinate of users
    D_position{n_cell}(:,2)=radius_d.*sin(phase)+BS_position{n_cell}(2); % y coodinate of users
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
        h_large{n_cell}(:,n_cell_1)=10.^(-(pl + shadow_factor)/10); % large scale fading coefficient
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

%% solve the maximum QoE and optimal k under every possible SNR for each user group
k=cell(N_cell,1); % the optimal k of each user group 
QoE_k=cell(N_cell,1); % the corresponding QoE
for n_cell=1:1:N_cell
    N_S=N(n_cell,1); % number of single-modal users in cell i
    N_Bi=N(n_cell,2); % number of bimodal users in cell i
    k{n_cell}=cell(N_Bi/2+N_S,1); % N_Bi/2+N_S user groups
    QoE_k{n_cell}=cell(N_Bi/2+N_S,1); % N_Bi/2+N_S user groups
    for q=1:1:N_Bi/2 % for each bimodal user group, solve the optimal k and corresponding QoE for every possible SNR combinations
        [k{n_cell}{q},QoE_k{n_cell}{q}]=Bi_SNR_k_QoE(SINR_Bi,VQA_table,H_Bi_text,H_Bi_image,K_Bi_text,K_Bi_image,bandwidth,para_Bi_text{n_cell}(q,:),para_Bi_image{n_cell}(q,:),G_th, w_phi{n_cell}(2*q-1:2*q));
    end
    for n_s=1:1:N_S % for each single-modal user, solve the optimal k and corresponding QoE for every possible SNR 
        n_ss=N_Bi/2+n_s; % the index of the single-modal user
        [k{n_cell}{n_ss},QoE_k{n_cell}{n_ss}]=Single_SNR_k_QoE(SINR_single,DeepSC_table,K_S,H_S,bandwidth,para_S{n_cell}(n_s,:),G_th,w_phi{n_cell}(N_Bi+n_s));
    end
end

%% all combinations of channels and power
com=cell(N_cell,2); % the first column: com_s; the second column: com_bi
for n_cell=1:1:N_cell
   N_D=N(n_cell,3); % number of all users in cell i
   channels=1:1:max(N_D,N_channels); % including real channels and virtual channels
   % all searching set of single modal users
   ind=fullfact([length(channels),length(P_range)]);
   com{n_cell,1}=[reshape(channels(ind(:,1)),1,length(ind(:,1)));reshape(P_range(ind(:,2)),1,length(ind(:,2)))];
   % all searching set of bimodal user pairs
   ind=fullfact([length(channels),length(channels),length(P_range),length(P_range)]);
   ind(find(ind(:,1)==ind(:,2)),:)=[]; % two channels can not be the same one
   com{n_cell,2}=[reshape(channels(ind(:,1)),1,length(ind(:,1)));reshape(channels(ind(:,2)),1,length(ind(:,2)));reshape(P_range(ind(:,3)),1,length(ind(:,3)));reshape(P_range(ind(:,4)),1,length(ind(:,4)))];
end

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

%% initial matching and states update
for n_cell=1:1:N_cell
    N_D=N(n_cell,3); % number of all users in cell i
    N_S=N(n_cell,1)+max(N_D,N_channels)-N_D; % number of single-modal users; add virtual users
    N_Bi=N(n_cell,2); % number of bimodal users in cell i
    ini_channel=randperm(max(N_D,N_channels)); % match each user a random channel
    for q=1:1:N_Bi/2 % for each bimodal user group, update the states of users and channels
        state_u{n_cell}{q}=[ini_channel(2*q-1:2*q);reshape(P_range(randi(length(P_range),2,1)),1,2);zeros(1,2);zeros(1,2)]; % 1st row: channels; 2nd row: powers; 3rd row: SINR; 4th row: QoE
        state_c{n_cell}(ini_channel(2*q-1:2*q),2)=[2*q-1,2*q];
    end
    for n_s=1:1:N_S % for each single-modal user, update the states of users and channels
        state_u{n_cell}{n_s+N_Bi/2}=[ini_channel(n_s+N_Bi);P_range(randi(length(P_range),1,1));0;0];
        state_c{n_cell}(ini_channel(n_s+N_Bi),2)=n_s+N_Bi; % update channel state
    end
end
% update QoE of user groups
for n_cell=1:1:N_cell
    N_D=N(n_cell,3); % number of all users in cell i
    N_S=N(n_cell,1)+max(N_D,N_channels)-N_D; % number of single-modal users; add virtual users
    N_Bi=N(n_cell,2); % number of bimodal users in cell i
    for q=1:1:N_Bi/2+N_S % all user groups updates their states, including the snr and QoE; update the corresponding channel states
        [state_u{n_cell}{q}(3,:),state_u{n_cell}{q}(4,:),state_c{n_cell}(state_u{n_cell}{q}(1,:),1)]=UpdateUC(state_u,q,n_cell,N_cell,state_c,QoE_k{n_cell},SINR_Bi,SINR_single,h,N,N_channels,P_noise);
    end 
end

%% iteration process
ite=0; % number of iterations
flag_swap=0; % record the number of swap operations
while 1 % until satisfy the end condition
    flag_ite=flag_swap; % save the last flag_swap
    ite=ite+1; % number of iterations
    for n_cell=1:1:N_cell % for each cell
        for q=1:1:N(n_cell,2)/2 % for bimodal user pairs
            for i_pre=1:1:length(com{n_cell,2}(1,:)) % search the resource set
                state_u_copy=state_u;
                state_c_copy=state_c;
               % swap channel pair: com{n_cell,2}(1,i_pre)--->state_u{n_cell}{q}(1,1)
               %                    com{n_cell,2}(2,i_pre)--->state_u{n_cell}{q}(1,2)
               %                    state_u{n_cell}{q}(1,1)--->index_d(com{n_cell,2}(1,i_pre))
               %                    state_u{n_cell}{q}(1,2)--->index_d(com{n_cell,2}(2,i_pre))
               % changed power level: com{n_cell,2}(3:4,i_pre) ---> state_u{n_cell}{q}(2,1:2)
               flag_y=1; % flag, to guarantee all related QoE is no less than the previous one 
               flag_h=0; % flag, to guarantee that exist one related QoE is larger than the previous one 
               
               state_u_copy{n_cell}{q}(2,:)=com{n_cell,2}(3:4,i_pre)'; % change power pair
               
               % try to swap the first channel in i_pre-th resource
               index_d1=state_c{n_cell}(com{n_cell,2}(1,i_pre),2); % the index of the user with the channel (com{n_cell,2}(1,i_pre))
               index_q=state_c{n_cell}(state_u_copy{n_cell}{q}(1,1),2); % the index of the user with the channel (state_u_copy{n_cell}{q}(1,1))
               % update the channel states
               state_c_copy{n_cell}(state_u_copy{n_cell}{q}(1,1),2)=index_d1; % change the user of state_u_copy{n_cell}{q}(1,1) to index_d1
               state_c_copy{n_cell}(com{n_cell,2}(1,i_pre),2)=index_q; % change the user of com{n_cell,2}(1,i_pre) to index_q
               if index_d1<=N(n_cell,2) % the user is a bimodal user
                   index_d1=ceil(index_d1/2); % the index of the bimodal user group
                   state_u_copy{n_cell}{index_d1}(1,(find(state_u{n_cell}{index_d1}(1,:)==com{n_cell,2}(1,i_pre))))=state_u{n_cell}{q}(1,1); % update the channel of index_d1
               else
                   index_d1=index_d1-N(n_cell,2)/2; % the user is a single-modal user; get the index of the user group
                   state_u_copy{n_cell}{index_d1}(1)=state_u{n_cell}{q}(1,1); % update the channel of index_d1
               end
               state_u_copy{n_cell}{q}(1,(find(state_u{n_cell}{q}(1,:)==state_u{n_cell}{q}(1,1))))=com{n_cell,2}(1,i_pre); % update the channel of q
               
               % copy state_u to avoid data mess
               state_u_temp=state_u_copy;
               
               % similar to index_d2; try to swap the second channel in i_pre-th resource
               index_d2=state_c_copy{n_cell}(com{n_cell,2}(2,i_pre),2);
               index_q=state_c_copy{n_cell}(state_u_copy{n_cell}{q}(1,2),2);
               % update the channel states
               state_c_copy{n_cell}(state_u_copy{n_cell}{q}(1,2),2)=index_d2;
               state_c_copy{n_cell}(com{n_cell,2}(2,i_pre),2)=index_q;
               if index_d2<=N(n_cell,2) % the user is a bimodal user
                   index_d2=ceil(index_d2/2);
                   state_u_copy{n_cell}{index_d2}(1,(find(state_u_temp{n_cell}{index_d2}(1,:)==com{n_cell,2}(2,i_pre))))=state_u_temp{n_cell}{q}(1,2);
               else
                   index_d2=index_d2-N(n_cell,2)/2;
                   state_u_copy{n_cell}{index_d2}(1)=state_u_temp{n_cell}{q}(1,2);
               end
               state_u_copy{n_cell}{q}(1,(find(state_u_temp{n_cell}{q}(1,:)==state_u_temp{n_cell}{q}(1,2))))=com{n_cell,2}(2,i_pre);
               
               % find the union of related user indexes in the current cell
               index=unique([q,index_d1,index_d2]); % find the unique user indexes
               for i=1:1:length(index) % for each related user
                   [state_u_copy{n_cell}{index(i)}(3,:),state_u_copy{n_cell}{index(i)}(4,:),state_c_copy{n_cell}(state_u_copy{n_cell}{index(i)}(1,:),1)]=UpdateUC(state_u_copy,index(i),n_cell,N_cell,state_c_copy,QoE_k{n_cell},SINR_Bi,SINR_single,h,N,N_channels,P_noise); % update the state of this user
                   if N(n_cell,3)<=N_channels % if the number of users is less than N_channel, we need to guarantee that the QoE of users do not decrease
                       flag_y=flag_y&&(state_u_copy{n_cell}{index(i)}(4,1)>=state_u{n_cell}{index(i)}(4,1)); % guarantee the QoE of this user is no less than the QoE before swapping
                       if flag_y==0 % if not satisify this condition, break out of this loop. Try the next possible resource. 
                          break
                       end
                       flag_h=flag_h || (state_u_copy{n_cell}{index(i)}(4,1)>state_u{n_cell}{index(i)}(4,1));
                   else % otherwise, guarantee the utility of channel will not decrease
                       for i_c=1:1:length(state_u_copy{n_cell}{index(i)}(1,:)) % for each channel of user i
                           flag_y=flag_y && (state_c_copy{n_cell}(state_u_copy{n_cell}{index(i)}(1,i_c),1)>=state_c{n_cell}(state_u_copy{n_cell}{index(i)}(1,i_c),1));% guarantee the QoE of this channel is no less than the QoE before swapping
                           if flag_y==0
                              break
                           end
                           flag_h=flag_h || (state_c_copy{n_cell}(state_u_copy{n_cell}{index(i)}(1,i_c),1)>state_c{n_cell}(state_u_copy{n_cell}{index(i)}(1,i_c),1));
                       end
                       if flag_y==0
                          break
                       end
                   end
               end
               % find the related user indexes in other cells (sharing the same channels)
               if flag_y==0
                  continue
               else
                   index_c=unique([state_u{n_cell}{q}(1,:),com{n_cell,2}(1:2,i_pre).']); % find the unique channels in this swap
                   index_c(index_c>N_channels)=[]; % delete virtual channels
                   for nn_cell=1:1:N_cell % for each cell
                       index=[]; % save the indexes of related users
                       if nn_cell~=n_cell % the cell should not be the current cell
                           for i_c=1:1:length(index_c) % for each channel
                               index_d=state_c{nn_cell}(index_c(i_c),2); % find the user index of this channel
                               if index_d<=N(nn_cell,2) % the user is a bimodal user
                                   index_d=ceil(index_d/2); % the index of this user group
                               else % the user is a single-modal user
                                   index_d=index_d-N(nn_cell,2)/2; % the index of this user group
                               end
                               index=[index,index_d]; % save the related users in nn_cell
                           end
                           for i=1:1:length(index) % solve the QoE of each related user
                               [state_u_copy{nn_cell}{index(i)}(3,:),state_u_copy{nn_cell}{index(i)}(4,:),state_c_copy{nn_cell}(state_u_copy{nn_cell}{index(i)}(1,:),1)]=UpdateUC(state_u_copy,index(i),nn_cell,N_cell,state_c_copy,QoE_k{nn_cell},SINR_Bi,SINR_single,h,N,N_channels,P_noise);
                               % decide whether the QoE of this user will not decrease; the code can be optimized
                               flag_y=flag_y&&(state_u_copy{nn_cell}{index(i)}(4,1)>=state_u{nn_cell}{index(i)}(4,1));
                               if flag_y==0
                                  break
                               end
                               flag_h=flag_h || (state_u_copy{nn_cell}{index(i)}(4,1)>state_u{nn_cell}{index(i)}(4,1));
                           end
                       end
                       if flag_y==0
                          break
                       end
                   end
               end
               if flag_y==1 && flag_h ==1 % the utility will increase, so the swap can be performed.
                   flag_swap=flag_swap+1;
                   state_u=state_u_copy;
                   state_c=state_c_copy;
                   break
               end
            end
        end
        
        % for single modal user
        for n_s=N(n_cell,2)/2+1:1:N(n_cell,2)/2+N(n_cell,1) % the index of the single-modal user
            for i_pre=1:1:length(com{n_cell,1}(1,:)) % searching the resource set of the single-modal task
                state_u_copy=state_u;
                state_c_copy=state_c;
               % swap channel pair: com{n_cell,1}(1,i_pre)--->state_u{n_cell}{q}(1)
               %                    state_u{n_cell}{q}(1)--->index_d(com{n_cell,1}(1,i_pre))
               % changed power level: com{n_cell,2}(2,i_pre) ---> state_u{n_cell}{q}(2)
               flag_y=1; % all related QoE is no less than the previous one 
               flag_h=0; % exist one related QoE larger than the previous one 
               
               state_u_copy{n_cell}{n_s}(2)=com{n_cell,1}(2,i_pre); % change the power pair
               
               % try to swap the channel in i_pre-th resource
               index_d1=state_c{n_cell}(com{n_cell,1}(1,i_pre),2);
               index_q=state_c{n_cell}(state_u_copy{n_cell}{n_s}(1),2);
               % update the channel states
               state_c_copy{n_cell}(state_u_copy{n_cell}{n_s}(1),2)=index_d1;
               state_c_copy{n_cell}(com{n_cell,1}(1,i_pre),2)=index_q;
               if index_d1<=N(n_cell,2) % the user is a bimodal user, and update the channel of index_d1
                   index_d1=ceil(index_d1/2);
                   state_u_copy{n_cell}{index_d1}(1,(find(state_u{n_cell}{index_d1}(1,:)==com{n_cell,1}(1,i_pre))))=state_u{n_cell}{n_s}(1);
               else % the user is a single-modal user, and then update the channel of index_d1
                   index_d1=index_d1-N(n_cell,2)/2;
                   state_u_copy{n_cell}{index_d1}(1)=state_u{n_cell}{n_s}(1);
               end
               state_u_copy{n_cell}{n_s}(1)=com{n_cell,1}(1,i_pre); % update the channel of n_s
               
               % find the union of related user indexes in the current cell
               index=unique([n_s,index_d1]);
               for i=1:1:length(index)
                   [state_u_copy{n_cell}{index(i)}(3,:),state_u_copy{n_cell}{index(i)}(4,:),state_c_copy{n_cell}(state_u_copy{n_cell}{index(i)}(1,:),1)]=UpdateUC(state_u_copy,index(i),n_cell,N_cell,state_c_copy,QoE_k{n_cell},SINR_Bi,SINR_single,h,N,N_channels,P_noise);
                   if N(n_cell,3)<=N_channels % QoE of users can not decrease
                       flag_y=flag_y&&(state_u_copy{n_cell}{index(i)}(4,1)>=state_u{n_cell}{index(i)}(4,1));
                       if flag_y==0
                          break
                       end
                       flag_h=flag_h || (state_u_copy{n_cell}{index(i)}(4,1)>state_u{n_cell}{index(i)}(4,1));
                   else % QoE of channels can not decrease
                       for i_c=1:1:length(state_u_copy{n_cell}{index(i)}(1,:))
                           flag_y=flag_y && (state_c_copy{n_cell}(state_u_copy{n_cell}{index(i)}(1,i_c),1)>=state_c{n_cell}(state_u_copy{n_cell}{index(i)}(1,i_c),1));
                           if flag_y==0
                              break
                           end
                           flag_h=flag_h || (state_c_copy{n_cell}(state_u_copy{n_cell}{index(i)}(1,i_c),1)>state_c{n_cell}(state_u_copy{n_cell}{index(i)}(1,i_c),1));
                       end
                   end
                   if flag_y==0
                       break
                   end
               end
               
               if flag_y==0
                   continue
               else
               % find the related user indexes in other cells (the channels)
                   index_c=unique([state_u{n_cell}{n_s}(1,:),com{n_cell,1}(1,i_pre)]);
                   index_c(index_c>N_channels)=[]; % delete virtual channels
                   for nn_cell=1:1:N_cell
                       index=[]; % save the indexes of related users
                       if nn_cell~=n_cell
                           for i_c=1:1:length(index_c)
                               index_d=state_c{nn_cell}(index_c(i_c),2);
                               if index_d<=N(nn_cell,2) % the user is a bimodal user
                                   index_d=ceil(index_d/2);
                               else
                                   index_d=index_d-N(nn_cell,2)/2;
                               end
                               index=[index,index_d]; 
                           end
                           for i=1:1:length(index) % solve the QoE of each related user
                               [state_u_copy{nn_cell}{index(i)}(3,:),state_u_copy{nn_cell}{index(i)}(4,:),state_c_copy{nn_cell}(state_u_copy{nn_cell}{index(i)}(1,:),1)]=UpdateUC(state_u_copy,index(i),nn_cell,N_cell,state_c_copy,QoE_k{nn_cell},SINR_Bi,SINR_single,h,N,N_channels,P_noise);
                               flag_y=flag_y&&(state_u_copy{nn_cell}{index(i)}(4,1)>=state_u{nn_cell}{index(i)}(4,1));
                               if flag_y==0
                                  break
                               end
                               flag_h=flag_h || (state_u_copy{nn_cell}{index(i)}(4,1)>state_u{nn_cell}{index(i)}(4,1));
                           end
                       end
                       if flag_y==0
                          break
                       end
                   end
               end
               if flag_y==1 && flag_h ==1 % exist one of QoE increasing
                   flag_swap=flag_swap+1;
                   state_u=state_u_copy;
                   state_c=state_c_copy;
                   break
               end
            end
        end
    end
    if flag_ite==flag_swap  % no swap operation in this time, end loop
        flag_ite
       break 
    end
end

%% save the results
for n_cell=1:1:N_cell
   for n_c=1:1:N_channels
       QoE_sum=QoE_sum+state_c{n_cell}(n_c,1);
   end
end

%% calculate the sum power consumption
% power_result=cell(N_cell,1);
% for n_cell=1:1:N_cell
%     N_S=N(n_cell,1);
%     N_Bi=N(n_cell,2);
%     power_result{n_cell}=zeros(1,N_S+N_Bi);
%     for q=1:1:N_Bi/2
%         if state_u{n_cell}{q}(4,1)==0
%            power_result{n_cell}(2*q-1:2*q)=[0,0];
%            continue
%         end
%        power_result{n_cell}(2*q-1:2*q)=state_u{n_cell}{q}(2,:);
%     end
%     for n_s=N_Bi/2+1:1:N_Bi/2+N_S
%        if state_u{n_cell}{n_s}(4,1)==0
%            power_result{n_cell}(n_s+N_Bi/2)=0;
%            continue
%        end
%        power_result{n_cell}(n_s+N_Bi/2)=state_u{n_cell}{n_s}(2);
%     end
% end
% 
% for n_cell=1:1:N_cell
%       P_sum=P_sum+sum(power_result{n_cell}); 
% end

% toc
end