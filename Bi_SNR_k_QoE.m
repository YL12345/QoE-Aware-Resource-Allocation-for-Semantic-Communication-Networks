function [k_Bi,QoE_Bi]=Bi_SNR_k_QoE(SINR_Bi,VQA_table,H_Bi_text,H_Bi_image,K_Bi_text,K_Bi_image,bandwidth,para_Bi_text,para_Bi_image,G_th, w_phi)
% solve the maximum QoE by optimizing k for this bimodal user pair, according to SNR pair
% Input: SINR_Bi: all possible combinations of SINR of two users
%       VQA_table: VQA_table=cell(length(K_Bi_image),length(K_Bi_text)); in each cell, row:snr of text user; column:snr of image user
%       H_Bi_text/H_Bi_image: semantic entropy of bimodal user
%       K_Bi_text/K_Bi_image: possible values of of number of semantic symbols for VQA
%       bandwidth
%       para_Bi_text/para_Bi_image: parameters related to QoE, including beta;phi;lamda;si
%       G_th=0.5; %the minimum score of phi and si for all users
%       w_phi: 1*2, the weights of phi for this bimodal user pair
% Output: SNR_Bi: 2*(length(snr_range)^2), all possible permutations of SNRs for Bimodal users
%         k_Bi: 2*(length(snr_range)^2),the corresponding optimal k pair for each SNR pair
%         QoE_Bi:  1*(length(snr_range)^2), the maximum QoE under each SNR pair

k_Bi=zeros(2,length(SINR_Bi));
QoE_Bi=zeros(1,length(SINR_Bi));

sinr_index=SINR_Bi/5+3; % obtain the SNR index, so as to find the semantic accuracy according to VQA_table

for i=1:1:length(sinr_index(1,:)) % for each snr pair
    sinr_t=sinr_index(1,i); % the index of sinr for text user
    sinr_i=sinr_index(2,i); % the index of sinr for image user
    QoE_q=zeros(length(K_Bi_text),length(K_Bi_image)); % save the QoE
    for k_t=1:1:length(K_Bi_text) % for each k_Bi_text
        for k_i=1:1:length(K_Bi_image) % for each k_Bi_image
            si=VQA_table{k_i,k_t}(sinr_t,sinr_i); % obtain the semantic similarity
            phi_t=H_Bi_text/(K_Bi_text(k_t)/bandwidth); % solve the phi for text user
            phi_i=H_Bi_image/(K_Bi_image(k_i)/bandwidth); % solve the phi for image user

            G_phi_t=1/(1+exp(para_Bi_text(1)*(para_Bi_text(2)-phi_t/1000))); %the score of phi for text user
            G_si_t=1/(1+exp(para_Bi_text(3)*(para_Bi_text(4)-si))); %the score of si for text user
            QoE_t=w_phi(1)*G_phi_t+(1-w_phi(1))*G_si_t; % QoE for text user

            G_phi_i=1/(1+exp(para_Bi_image(1)*(para_Bi_image(2)-phi_i/1000))); % the score of phi for image user
            G_si_i=1/(1+exp(para_Bi_image(3)*(para_Bi_image(4)-si)));% the score of si for image user
            QoE_i=w_phi(2)*G_phi_i+(1-w_phi(2))*G_si_i; % QoE for image user
            if G_phi_t<G_th || G_si_t<G_th || G_phi_i<G_th || G_si_i<G_th % if the scores is less than the threshold, set the QoE to 0
                QoE_q(k_t,k_i)=0;
            else % otherwise, the QoE is the sum of QoE_t and QoE_i
                QoE_q(k_t,k_i)=QoE_t+QoE_i;
            end
        end
    end
    [temp, index_t]=max(QoE_q); %find the maximum QoE for each k_i
    [temp, index_i]=max(temp); %find the maximum QoE
    QoE_Bi(i)=temp; %save the maximum QoE of all (k_t,k_i)
    k_Bi(:,i)=[index_t(index_i),index_i]; %save the index of the corresponding k pair
end


end