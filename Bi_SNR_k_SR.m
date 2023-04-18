function [k_Bi,SR_Bi]=Bi_SNR_k_SR(SINR_Bi,VQA_table,H_Bi_text,H_Bi_image,K_Bi_text,K_Bi_image,bandwidth,SR_th_Bit,SR_th_Bii)
% solve the maximum Semantic transmission rate by optimizing k for this bimodal user pair, according to SNR pair
% Input: SINR_Bi: all possible combinations of SINR of two users
%       VQA_table: VQA_table=cell(length(K_Bi_image),length(K_Bi_text)); in each cell, row:snr of text user; column:snr of image user
%       H_Bi_text/H_Bi_image: semantic entropy of bimodal user
%       K_Bi_text/K_Bi_image: possible values of of number of semantic symbols for VQA
%       bandwidth
%       SR_th_Bit/SR_th_Bii; %the threshold of SR for all users
% Output: SNR_Bi: 2*(length(snr_range)^2), all possible permutations of SNRs for Bimodal users
%         k_Bi: 2*(length(snr_range)^2),the corresponding optimal k pair for each SNR pair
%         QoE_Bi:  1*(length(snr_range)^2), the maximum QoE under each SNR pair

k_Bi=zeros(2,length(SINR_Bi));
SR_Bi=zeros(1,length(SINR_Bi));

sinr_index=SINR_Bi/5+3; %so as to find the semantic accuracy according to sem_table

for i=1:1:length(sinr_index(1,:)) % for each sinr pair
    sinr_t=sinr_index(1,i); % the index of sinr for text user
    sinr_i=sinr_index(2,i); % the index of sinr for image user
    SR_q=zeros(length(K_Bi_text),length(K_Bi_image)); % save the S-R of different k pair
    si_temp=zeros(length(K_Bi_text),length(K_Bi_image)); % save si
    % find the optimal k using exhausted searching
    for k_t=1:1:length(K_Bi_text) % for each k_t
        for k_i=1:1:length(K_Bi_image) % for each k_i
            si=VQA_table{k_i,k_t}(sinr_t,sinr_i); % obtain si
            si_temp(k_t,k_i)=si;
            SR_t=H_Bi_text*si/(K_Bi_text(k_t)/bandwidth); % solve S-R for text transmission user
            SR_i=H_Bi_image*si/(K_Bi_image(k_i)/bandwidth); % solve S-R for image transmission user
            % decide whether the S-R is larger than the threshold; if no, set it to 0
            if SR_t<SR_th_Bit(2)*SR_th_Bit(1) || SR_i<SR_th_Bii(2)*SR_th_Bii(1) || si<SR_th_Bit(2) || si<SR_th_Bii(2)
                SR_q(k_t,k_i)=0;
            else % if yes, the S-R of this bimodal user group is the sum of SR_t and SR_i
                SR_q(k_t,k_i)=SR_t+SR_i;
            end
        end
    end
    [temp, index_t]=max(SR_q); %find the maximum SR for each k_i
    [temp, index_i]=max(temp); %find the maximum sum SR
    SR_Bi(i)=temp; %save the maximum SR of all (k_t,k_i)
    k_Bi(:,i)=[index_t(index_i),index_i]; %save the index of the corresponding k pair
end


end