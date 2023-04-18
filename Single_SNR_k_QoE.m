function [k_single,QoE_single]=Single_SNR_k_QoE(SINR_single,DeepSC_table,K_S,H_S,bandwidth,para_S,G_th,w_phi)
% solve the optimal k and maximum QoE under different snr for each single-modal user
% Input: SNR_single: the range of SNR for single-modal users
%        DeepSC_table: the performance table of single-modal users; row: number of symbols; column:snr
%        K_S: the range of k for single-modal users
%        H_S: semantic entropy of single-modal user
%        bandwidth: the bandwidth of each channel
%       para_S: parameters related to QoE, including beta_S;phi_S;lamda_S;si_S
%       G_th=0.5; %the minimum score of phi and si for all users
%       w_phi: the weights of phi for this single-modal user
% Output: k_single: the corresponding optimal k for each SNR
%         QoE_single: the maximum QoE under each SNR

k_single=zeros(1,length(SINR_single)); % save the results
QoE_single=zeros(1,length(SINR_single));

sem_table = DeepSC_table; % data for DeepSC

sinr_index=SINR_single-min(SINR_single)+1; %the snr index, so as to find the semantic accuracy according to sem_table

for i=1:1:length(SINR_single) % for each possible snr
    QoE_temp=zeros(1,length(K_S));
   for k=1:1:length(K_S) % for each k
       si=sem_table(k,sinr_index(i)); %the semantic accuracy for different k
       phi=H_S/(k/bandwidth); %semantic emit rate for every k
       G_phi=1/(1+exp(para_S(1)*(para_S(2)-phi/1000))); %the score of phi; phi/1000 is to change the unit to ksut/s
       G_si=1/(1+exp(para_S(3)*(para_S(4)-si))); %the score of si
       QoE=w_phi*G_phi+(1-w_phi)*G_si; %QoE
       if G_si<G_th || G_phi<G_th % if the constrants can not be satisfied, set QoE and k as 0
          QoE_temp(k)=0;
       else
           QoE_temp(k)=QoE;
       end
   end
   [QoE, index]=max(QoE_temp); %find the maximum QoE
   QoE_single(i)=QoE;
   k_single(i)=K_S(index);
end


end