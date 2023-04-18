function [k_single,SR_single]=Single_SNR_k_SR(SINR_single,DeepSC_table,K_S,H_S,bandwidth,SR_th_S)
% for each single-modal user, solve the maximum SR by optimizing k, according to SNR
% Input: SNR_single: the range of SNR for single-modal users
%        DeepSC_table: the performance table of single-modal users; row: number of symbols; column:snr
%        K_S: the range of k for single-modal users
%        H_S: semantic entropy of single-modal user
%        bandwidth: the bandwidth of each channel
%       SR_th_S; 1*2, the minimum SR for all users
% Output: k_single: the corresponding optimal k for each SNR
%         SR_single: the maximum QoE under each SNR

k_single=zeros(1,length(SINR_single));
SR_single=zeros(1,length(SINR_single));

% the performance table for DeepSC
sem_table = DeepSC_table;

sinr_index=SINR_single-min(SINR_single)+1; %so as to find the semantic accuracy according to sem_table

for i=1:1:length(SINR_single) % for each sinr
    SR_temp=zeros(1,length(K_S));
   for k=1:1:length(K_S) % for each k; using exhausted searching to find the optimal results
       si=sem_table(k,sinr_index(i)); %the semantic accuracy for different k
       SR=H_S*si/(k/bandwidth); %semantic transmission rate for every k

       if si<SR_th_S(2) || SR<SR_th_S(2)*SR_th_S(1) % if the constrants can not be satisfied, set SRand k to 0
          SR_temp(k)=0;
       else
           SR_temp(k)=SR;
       end
   end
   [SR, index]=max(SR_temp); %find the maximum SR
   SR_single(i)=SR;
   k_single(i)=K_S(index);
end


end