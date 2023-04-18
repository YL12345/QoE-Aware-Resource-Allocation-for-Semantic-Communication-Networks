function [SINR,QoE,QoE_c]=UpdateUC(state_u,index_ug,n_cell,N_cell,state_c,QoE_k,SINR_Bi,SINR_single,h,N,N_channels,P_noise)
% update SINR and QoE of user group (index_ug) and QoE of corresponding channels
% Input: state_u: cell: N_cell*1; cell: (N_Bi/2+N_S)*1; 4*2/4*1: 1st row: channels; 2nd row: powers; 3rd row: SINR; 4th row: QoE
%        index_ug: the index of the current user group. bimodal if <=N_Bi/2; otherwise,singlemodal user 
%        n_cell: the current cell
%        N_cell: the number of cells
%        state_c: cell: N_cell*1; max(N_D,N_channels)*2; row: channels; 1st column: the current QoE; 2rd column: index of user groupstates of channels, used for finding the user indexes with the same channel in other cells
%        QoE_k: cell:(N_Bi/2+N_S)*1; 1*length(SINR_Bi/SINR_single)
%        SINR_Bi: 2*length(snr_range)^2; all possible combinations of SINR of two users
%        SINR_single: 1*length(snr_range); all possible SINR of single-modal users
%        h: cell: N_cell*N_channels; Each cell: N_cell*Nr*N_D; including large scale fading and small scale fading
%        N: N_cell*3; number of users; the first column: N_S; the second column: N_Bi; the third column: N_D=N_S+N_Bi
%        N_channels: number of channels
%        P_noise: the noise power
% Output: SINR: 1*2/1*1, the SINR of index_ug
%         QoE: 1*2/1*1, the QoE of index_ug
%         QoE_c: 2*1/1*1,the QoE of channels

if length(state_u{n_cell}{index_ug}(1,:))==2 % bimodal user pair
    if state_u{n_cell}{index_ug}(1,1)>N_channels || state_u{n_cell}{index_ug}(1,2)>N_channels  % it is a virtual channel. Set the output to 0.
        SINR=zeros(1,2);
        QoE=zeros(1,2);
        QoE_c=zeros(2,1);
        return
    end
    n_d1=2*index_ug-1; % the index of the user 1
    channel_1=state_u{n_cell}{index_ug}(1,1); % the corresponding channel index
    I1=0; % the interference of n_d1
    n_d2=2*index_ug;% the index of the user 2
    channel_2=state_u{n_cell}{index_ug}(1,2);% the corresponding channel index
    I2=0;% the interference of n_d2
    H1=h{n_cell,channel_1}(n_cell,:,n_d1).'; % the channel coefficient matrix for MRC detection of n_d1
    H2=h{n_cell,channel_2}(n_cell,:,n_d2).'; % the channel coefficient matrix for MRC detection of n_d2
    %calculate the interference experienced at each user
    for nn_cell =1:1:N_cell % for each cell
       if nn_cell ~=n_cell % if it is not n_cell
          ind1=state_c{nn_cell}(channel_1,2); % find the index of interference user in nn_cell over channel_1
          if ind1<=N(nn_cell,3) % it is not a virtual user
              H_i=h{nn_cell,channel_1}(n_cell,:,ind1).'; % get H_i of the interference user
              if ind1>N(nn_cell,2) % single modal user
                  ind1_ug=ind1-N(nn_cell,2)/2; % get the index
                  I1=I1+state_u{nn_cell}{ind1_ug}(2)*(abs(H1'*H_i))^2; % calculate the interference
              else % bimodal user
                  ind1_ug=ceil(ind1/2); % the index of bimodal user pair
                  I1=I1+state_u{nn_cell}{ind1_ug}(2,ind1-2*(ind1_ug-1))*(abs(H1'*H_i))^2; % solve the interference
              end
          end
          ind2=state_c{nn_cell}(channel_2,2); %find the index of interference user in nn_cell over channel_2
          if ind2<=N(nn_cell,3) % it is not a virtual user
              H_i=h{nn_cell,channel_2}(n_cell,:,ind2).';
              if ind2>N(nn_cell,2) % single modal user
                  ind2_ug=ind2-N(nn_cell,2)/2;
                  I2=I2+state_u{nn_cell}{ind2_ug}(2)*(abs(H2'*H_i))^2;
              else
                  ind2_ug=ceil(ind2/2); % the index of bimodal user pair
                  I2=I2+state_u{nn_cell}{ind2_ug}(2,ind2-2*(ind2_ug-1))*(abs(H2'*H_i))^2;
              end
          end
       end
    end
    % SINR 
    SINR1=(state_u{n_cell}{index_ug}(2,1)*(H1'*H1)^2)/((norm(H1',2)^2)*P_noise+I1); % solve the SINR
    sinr_db1=round(10*log10(SINR1)); % (dB)
    SINR2=(state_u{n_cell}{index_ug}(2,2)*(H2'*H2)^2)/((norm(H2',2)^2)*P_noise+I2);
    sinr_db2=round(10*log10(SINR2));
    SINR=[sinr_db1,sinr_db2];
    sinr_temp=SINR.'-SINR_Bi;
    temp=max(find(sinr_temp(1,:)>=0&sinr_temp(2,:)>=0)); % find the index of SINR in SINR_Bi
    if isempty(temp) % no such a sinr pair in SINR_Bi, which means the sinr is to small to transmit data. Set the result to 0.
       QoE=zeros(1,2);
       QoE_c=zeros(2,1);
       return
    else
       QoE=QoE_k{index_ug}(temp)/2*ones(1,2); % find the QoE
       QoE_c=QoE'; % save the QoE to the states of the correponding channels
    end
    
else % single modal user
    if state_u{n_cell}{index_ug}(1)>N_channels || index_ug>N(n_cell,2)/2+N(n_cell,1)  % it is a virtual channel or a virtual user
        SINR=0;
        QoE=0;
        QoE_c=0;
        return
    end
    n_d=index_ug+N(n_cell,2)/2; % the index of the user
    channel=state_u{n_cell}{index_ug}(1); % the corresponding channel
    I=0; % the interference
    H=h{n_cell,channel}(n_cell,:,n_d).'; % the channel coefficient matrix for MRC detection
    %calculate the interference experienced at each user
    for nn_cell =1:1:N_cell 
       if nn_cell ~=n_cell
          ind=state_c{nn_cell}(channel,2); %find the index of interference user in nn_cell over channel_1
          if ind<=N(nn_cell,3) % it is not a virtual user
              H_i=h{nn_cell,channel}(n_cell,:,ind).';
              if ind>N(nn_cell,2) % single modal user
                  ind_ug=ind-N(nn_cell,2)/2;
                  I=I+state_u{nn_cell}{ind_ug}(2)*(abs(H'*H_i))^2;
              else
                  ind_ug=ceil(ind/2); % the index of bimodal user pair
                  I=I+state_u{nn_cell}{ind_ug}(2,ind-2*(ind_ug-1))*(abs(H'*H_i))^2;
              end
          end
       end
    end
    % SINR 
    SINR=(state_u{n_cell}{index_ug}(2,1)*(H'*H)^2)/((norm(H',2)^2)*P_noise+I);
    SINR=round(10*log10(SINR));
    sinr_temp=SINR-SINR_single;
    temp=max(find(sinr_temp>=0)); % the index of SINR in SINR_Bi
    if isempty(temp)
        QoE=0;
        QoE_c=0;
    else
        QoE=QoE_k{index_ug}(temp);
        QoE_c=QoE;
    end
end

end