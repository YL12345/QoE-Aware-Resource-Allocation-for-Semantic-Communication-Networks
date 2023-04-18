clc
clear
tic
% main function
% 3 cells; uplink transmission; 
% semantic-aware networks
% two applications: DeepSC for text transmission and VQA task; single modal task and bimodal task

%% cell parameters
N_cell = 3; 
ISD=1000; % the distance between two BSs
radius=ISD/sqrt(3); % cell radius

% the position of each BS
BS_position=cell(N_cell,1); 
BS_position{1}=[0,0]; %the position of BS1:(0,0)
BS_position{2}=[ISD,0]; % the position of BS2:(1000,0)
BS_position{3}=[ISD/2,ISD*sqrt(3)/2];
% BS_position{4}=[ISD/2,-ISD*sqrt(3)/2];
% BS_position{5}=[-ISD,0];
% BS_position{6}=[-ISD/2,ISD*sqrt(3)/2];
% BS_position{7}=[-ISD/2,-ISD*sqrt(3)/2];

%% appilication parameters
% single-modal task
H_S=4; % semantic entropy of DeepSC
K_S=1:1:20; % possible values of the number of semantic symbols for DeepSC
load('DeepSC_table'); % row: number of symbols(1:1:20); column:snr(-10:1:20)

% bimodal task
H_Bi_text=2; % semantic entropy of text user of VQA
H_Bi_image=4*197; % semantic entropy of image user of VQA
K_Bi_text=[4,8,12,16,20]/2; % possible values of number of semantic symbols for text user in VQA; unit: symbols/word
K_Bi_image=[4,8,16,24,32]/2*197; % possible values of number of semantic symbols for image user in VQA; unit: symbols/image
snr_range=-10:5:20; % the range of snr for each bimodal user
ind=fullfact([length(snr_range),length(snr_range)]);
SINR_Bi=[snr_range(ind(:,1));snr_range(ind(:,2))]; % the snr combinations
load('VQA_table'); %VQA_table=cell(length(K_Bi_image),length(K_Bi_text)); in each cell, row:snr of text user; column:snr of image user

G_th=0.5; %the minimum score of phi and si for all users

%% channel parameters
N_channels=6; % the number of channels
shadow_factor=6; % shadowing factor(dB)
bandwidth=180000; %channel bandwidth is 180kHz

Nr=2;% the number of receive antennas equiped at each BS
Nt=1; % the number of transmit antennas equiped at each user

P_range_db=[-10,-5,0,5,10,15,20]; % the transmit power (dBm)
P_range=10.^((P_range_db)/10); % the transmit power(mW)
P_noise=180000 * 10^(-17.4); % noise power with 180kHz bandwidth (mW)(PSD:-174dBm/Hz)

I_th=0; % interference power threshold£¨not used in the work£©

%% monte-carlo simulation
% the method_index of the schemes:
% 1. exhausted searching; (not used in the work)
% 2. random matching--randomly match the channel and power for each user 
% 3. matching; randomly initialize the channels and power for each user (not used in the paper)
% 4. fix power(P_max) + PropoedMatching;(not used in the paper)
% 5. initial power is set as the smallest value to ensure higher energy efficiency; channel: random matching (Proposed algorithm)
% 6. single cell optimization; comparison with the scheme without considering inter-cell interference
% 7. the scheme where sum S-SR is optimized

% 2,5,6,7
method_index=5; % the index of simulated method

% channels_range=[5]; % number of channels
channels_range=1:1:7; % number of channels
% QoE_channel=zeros(1,length(channels_range)); % save the sum overall QoE result of x-axis
% P_channel=zeros(1,length(channels_range));
% up_channel=zeros(1,length(channels_range));

% G_th_range=0.3:0.1:0.9;
% QoE_channel=zeros(1,length(G_th_range));

N_index=1:1:3; % the setting for different numbers of users

Monte=1;
QoE_result_channel=[];

for i=1:1:Monte
    QoE_result=[]; % save the overall QoE in each time    
    P_result=[];
    up_result=[]; % save the upper bound
    for N_channels=channels_range
%     for n_index=N_index % different numbers of users
        % number of users in each cell
         N=zeros(N_cell,3); % first column: N_S; the second column: N_Bi,the third column: N_D=N_S+N_Bi

%         switch n_index % different n_index means different setting 
%             case 1
%                 N_S_sum=6; % the sum of single modal users in N_cell cells
%                 N_Bi_sum=6; % the sum of bimodal user pairs in N_cell cells
%             case 2
%                 N_S_sum=6; % the sum of single modal users in N_cell cells
%                 N_Bi_sum=9; % the sum of bimodal user pairs in N_cell cells
%             case 3
%                 N_S_sum=9; % the sum of single modal users in N_cell cells
%                 N_Bi_sum=9; % the sum of bimodal user pairs in N_cell cells
%         end

        N_S_sum=6; % default setting for the number of single-modal users; the sum of single modal users in N_cell cells
        N_Bi_sum=6; % default setting for the number of bimodal user groups; the sum of bimodal user pairs in N_cell cells
        cell_S_index=randi(N_cell,1,N_S_sum); % the number of single-modal users in each cell
        cell_Bi_index=randi(N_cell,1,N_Bi_sum);% the number of bimodal user groups in each cell
        
        % the users in every cell; meanwhile obtain the upper bound
        up_cell=0;
        for n_cell=1:1:N_cell % for each cell
            N_S=sum(cell_S_index==n_cell); % the number of single-modal users in n_cell
            N_Bi=sum(cell_Bi_index==n_cell)*2; % the number of bimodal users in n_cell
            N(n_cell,:)=[N_S,N_Bi,N_Bi+N_S]; % the users in n_cell
           if N_Bi+N_S>N_channels % if the number of users is larger than the number of channels, there are totally at most N_channels users can be served.
               up_cell=up_cell+N_channels;
           else % else there are totally N_Bi+N_S users can be served.
               up_cell=up_cell+N_Bi+N_S;
           end
        end
        up_result=[up_result,up_cell]; % save the number of served users in this time

        switch method_index
            case 1 % exhausted (not used in the paper)
                [QoE]=Exhausted(BS_position,SINR_single,DeepSC_table,SINR_Bi,VQA_table,N,radius,N_cell,N_channels,shadow_factor,Nr,P_range,P_noise,I_th,H_S,K_S,bandwidth,H_Bi_text,H_Bi_image,K_Bi_text,K_Bi_image,G_th);
            case 2 % random matching
                [QoE]=Random(BS_position,SINR_single,DeepSC_table,SINR_Bi,VQA_table,N,radius,N_cell,N_channels,shadow_factor,Nr,P_range,P_noise,I_th,H_S,K_S,bandwidth,H_Bi_text,H_Bi_image,K_Bi_text,K_Bi_image,G_th);
            case 3 % matching (not used in the paper)
                [QoE,state_u]=Matching(BS_position,SINR_single,DeepSC_table,SINR_Bi,VQA_table,N,radius,N_cell,N_channels,shadow_factor,Nr,P_range,P_noise,I_th,H_S,K_S,bandwidth,H_Bi_text,H_Bi_image,K_Bi_text,K_Bi_image,G_th);
            case 4 % FixP+matching (input P_range as only a value) (not used in the paper)
%                 P_range=10.^(20/10);
%                 P_range=10.^(0/10);
                [QoE,state_u,P]=MinPowerMatching(BS_position,SINR_single,DeepSC_table,SINR_Bi,VQA_table,N,radius,N_cell,N_channels,shadow_factor,Nr,P_range,P_noise,I_th,H_S,K_S,bandwidth,H_Bi_text,H_Bi_image,K_Bi_text,K_Bi_image,G_th);
            case 5 % ProposedAlgorithm; initial power is set as the minimum value to minimize the energy consumption
                [QoE,state_u,P]=MinPowerMatching(BS_position,SINR_single,DeepSC_table,SINR_Bi,VQA_table,N,radius,N_cell,N_channels,shadow_factor,Nr,P_range,P_noise,I_th,H_S,K_S,bandwidth,H_Bi_text,H_Bi_image,K_Bi_text,K_Bi_image,G_th);
            case 6 % single cell optimization; comparison with the scheme without considering inter-cell interference
                [QoE,state_u]=SCMatching(BS_position,SINR_single,DeepSC_table,SINR_Bi,VQA_table,N,radius,N_cell,N_channels,shadow_factor,Nr,P_range,P_noise,I_th,H_S,K_S,bandwidth,H_Bi_text,H_Bi_image,K_Bi_text,K_Bi_image,G_th);
            case 7 % the scheme where sum S-SR is optimized
                [QoE,state_u]=MaxSRMatching(BS_position,SINR_single,DeepSC_table,SINR_Bi,VQA_table,N,radius,N_cell,N_channels,shadow_factor,Nr,P_range,P_noise,I_th,H_S,K_S,bandwidth,H_Bi_text,H_Bi_image,K_Bi_text,K_Bi_image,G_th);
        end
        QoE_result=[QoE_result, QoE]; % save the results of QoE in this time
%         P_result=[P_result, P];
    end
    QoE_result_channel=[QoE_result_channel; QoE_result]; % save the QoE results when the number of channels is N_channels
%     up_channel=up_channel+up_result; 
%     P_channel=P_channel+P_result;
end

QoE_channel=sum(QoE_result_channel);
QoE_channel=QoE_channel/Monte;
figure;
plot(channels_range, QoE_channel);
% plot(N_index, QoE_channel);

% up_channel=up_channel/Monte;
% figure;
% plot(channels_range, up_channel);

% P_channel=P_channel/Monte;
% figure;
% plot(channels_range, P_channel);

toc



