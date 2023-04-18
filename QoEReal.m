function [QoE_real]=QoEReal(K_Bi_text, K_Bi_image, para_Bi_text, para_Bi_image, H_Bi_text, H_Bi_image, G_th, w_phi, para_S, bandwidth, H_S,N, channel_result, k_result, power_result, N_cell, P_noise, h)
% Calculate the actual QoE 
% Input: N_S: number of single-modal users
%        channel_result:N_cell*N_D; content:index of channels
%        k_result: N_cell*N_D; content: k of all users; 0 means no channel assigned to this user
%        power_result: N_cell*N_D; content: transmit power of all uses
%        N_cell: number of cells
%        P_noise: the noise power (mW)
%        h: cell: N_cell*N_channels; Each cell: N_cell*Nr*N_D; including large scale fading and small scale fading
% Output: QoE_real: N_cell*N_D; content QoE of all users

QoE_real = cell(N_cell, 1);

% the performance table of DeepSC for the text transmission task
snr_range_S=-10:1:20;
load('DeepSC_table'); % row: number of symbols; column:snr
sem_table = DeepSC_table;

% the performance table of DeepSC-VQA for the VQA task
snr_range_Bi=-10:5:20;
load('VQA_table'); %VQA_table=cell(length(K_Bi_image),length(K_Bi_text)); in each cell, row:snr of text user; column:snr of image user

for n_cell=1:1:N_cell % for each cell
    N_D=N(n_cell,3);
    N_S=N(n_cell,1);
    N_Bi=N(n_cell,2);
    sinr_db=zeros(1,N_D);
    for n_d=1:1:N_D % for each user
        if channel_result{n_cell}(n_d)~=0 % if this user is assigned a channel
           I=0; % the interference
           channel_i= channel_result{n_cell}(n_d); % the index of the channel
           H=h{n_cell,channel_i}(n_cell,:,n_d).'; % the channel coefficient matrix for MRC detection
           %calculate the interference experienced at each user
           for nn_cell =1:1:N_cell 
               if nn_cell ~=n_cell
                  in_d_i=find(channel_result{nn_cell}==channel_i); %find the index of interference user in nn_cell over channel_i
                  if isempty(in_d_i)
                      continue
                  else
                    H_i=h{nn_cell,channel_i}(n_cell,:,in_d_i).';
                    I=I+power_result{nn_cell}(in_d_i)*(abs(H'*H_i))^2;
                  end
               end
           end
           % SINR 
           SINR=(power_result{n_cell}(n_d)*(H'*H)^2)/((norm(H',2)^2)*P_noise+I);
           sinr_db(n_d)=round(10*log10(SINR));
        end
    end
    % QoE of single-modal users
    for nn_s=1:1:N_S 
        n_s=nn_s+N_Bi;
       if k_result{n_cell}(n_s)~=0 % if k is not 0
          if sinr_db(n_s) < min(snr_range_S)
               QoE_real{n_cell}(n_s)=0; % obtain the QoE
           else
               if sinr_db(n_s)>max(snr_range_S) % if SINR is larger than max(snr_range_S), set it to max(snr_range_S)
                   sinr_db(n_s)=max(snr_range_S);
               end
               sinr_index=sinr_db(n_s)-min(snr_range_S)+1;
               si=sem_table(k_result{n_cell}(n_s),sinr_index); %the semantic accuracy 
               phi=H_S/(k_result{n_cell}(n_s)/bandwidth); %semantic rate 
               G_phi=1/(1+exp(para_S{n_cell}(nn_s,1)*(para_S{n_cell}(nn_s,2)-phi/1000))); %the score of phi; phi/1000 is to change the unit to ksut/s
               G_si=1/(1+exp(para_S{n_cell}(nn_s,3)*(para_S{n_cell}(nn_s,4)-si))); %the score of si
               QoE=w_phi{n_cell}(n_s)*G_phi+(1-w_phi{n_cell}(n_s))*G_si; %QoE
               if G_si<G_th || G_phi<G_th % if the constrants can not be satisfied, QoE and k are 0
                  QoE_real{n_cell}(n_s)=0;
               else
                  QoE_real{n_cell}(n_s)=QoE;
               end                 
           end 
       end    
    end
    % QoE of bimodal users
    for q=1:1:N_Bi/2
       text_i=2*q-1; % the index of text transmission user
       image_i=2*q; % the index of image transmission user
       if k_result{n_cell}(image_i)~=0 && k_result{n_cell}(text_i)~=0
           if sinr_db(text_i) < min(snr_range_Bi) || sinr_db(image_i) < min(snr_range_Bi) % if the sinr is too small, QoE is 0
               QoE_real{n_cell}(text_i)=0;
               QoE_real{n_cell}(image_i)=0;
               continue
           else % if SINR is larger than the maximum, set it to the maximum value
               if sinr_db(text_i)>max(snr_range_Bi)
                   sinr_db(text_i)=max(snr_range_Bi);
               end
               if sinr_db(image_i)>max(snr_range_Bi)
                   sinr_db(image_i)=max(snr_range_Bi);
               end
           end
           sinr_index_t=floor(sinr_db(text_i)/5)+3; % the index of sinr of text transmission user
           sinr_index_i=floor(sinr_db(image_i)/5)+3;% the index of sinr of image transmission user
           si=VQA_table{k_result{n_cell}(image_i),k_result{n_cell}(text_i)}(sinr_index_t,sinr_index_i); % the semantic accuracy
           phi_t=H_Bi_text/(K_Bi_text(k_result{n_cell}(text_i))/bandwidth); % semanti rate of text transmission user
           phi_i=H_Bi_image/(K_Bi_image(k_result{n_cell}(image_i))/bandwidth); % semanti rate of image transmission user

           % the QoE of text transmission user
           G_phi_t=1/(1+exp(para_Bi_text{n_cell}(q,1)*(para_Bi_text{n_cell}(q,2)-phi_t/1000))); %the score of phi
           G_si_t=1/(1+exp(para_Bi_text{n_cell}(q,3)*(para_Bi_text{n_cell}(q,4)-si))); %the score of si
           QoE_t=w_phi{n_cell}(text_i)*G_phi_t+(1-w_phi{n_cell}(text_i))*G_si_t;
            % the QoE of image transmission user
           G_phi_i=1/(1+exp(para_Bi_image{n_cell}(q,1)*(para_Bi_image{n_cell}(q,2)-phi_i/1000)));
           G_si_i=1/(1+exp(para_Bi_image{n_cell}(q,3)*(para_Bi_image{n_cell}(q,4)-si)));
           QoE_i=w_phi{n_cell}(image_i)*G_phi_i+(1-w_phi{n_cell}(image_i))*G_si_i;
           if G_phi_t<G_th || G_si_t<G_th || G_phi_i<G_th || G_si_i<G_th % the scores need to be larger than the threshold
               QoE_real{n_cell}(text_i)=0;
               QoE_real{n_cell}(image_i)=0;
               continue
           else
               QoE_real{n_cell}(text_i)=QoE_t;
               QoE_real{n_cell}(image_i)=QoE_i;
           end 
       end
    end
end


end