function [ p_d_final,p_i_r_final,p_r_final,p_i_d_final ] = Channel_Analusis( X_Tx,Y_Tx,X_Rx,Y_Rx,X_A,Y_A,Number_Rx,MapLength,G_vec )
%CHANNAL 此处显示有关此函数的摘要
%   此处显示详细说明
Time_Measure = true;
if(Time_Measure)
    tic;
end

%% Parameters
sigma_shadowing=4;
B=30;
K=30;
gammaT_dB=10;
gammaT=10^(gammaT_dB/10);
noise_power_dB=-144;
noise_power=10^(noise_power_dB/10);

% Rx(CLPC)
X_Rx_CLPC = rand()*2*MapLength - MapLength;
Y_Rx_CLPC = rand()*2*MapLength - MapLength;

%% Main Process
DTR = zeros(Number_Rx,1);
DAR = zeros(Number_Rx,1);
for NoRx = 1:Number_Rx
    DTR(NoRx) = sqrt((X_Tx-X_Rx(NoRx))^2 + (Y_Tx-Y_Rx(NoRx))^2);
	DAR(NoRx) = sqrt((X_A-X_Rx(NoRx))^2 + (Y_A-Y_Rx(NoRx))^2);
end
DTA = sqrt((X_A-X_Tx)^2 + (Y_A-Y_Tx)^2);

gpTR_dB_vec = -128.1 -37.6*log10(DTR./(10^3));
gpTR_vec =(10.^(gpTR_dB_vec./10));
gpAR_dB_vec = -128.1 -37.6*log10(DAR./(10^3));
gpAR_vec =(10.^(gpAR_dB_vec./10));
gpTA_dB = -128.1 -37.6*log10(DTA/(10^3));
gpTA =(10^(gpTA_dB/10));

gsTA_dB = normrnd(0,sigma_shadowing);
gsTA = 10^(gsTA_dB/10);
gsAR_dB_vec = normrnd(0,sigma_shadowing,1,Number_Rx);
gsAR_vec=10.^(gsAR_dB_vec./10);
gsTR_dB_vec = normrnd(0,sigma_shadowing,1,Number_Rx);
gsTR_vec=10.^(gsTR_dB_vec./10);

h_TA_s=sqrt(1/2)*(randn(B,Number_Rx)+1i*randn(B,Number_Rx));
% hh_TA_s=sqrt(1/2)*(randn(B,Number_Rx)+1i*randn(B,Number_Rx));
h_TR_s=sqrt(1/2)*(randn(B,Number_Rx)+1i*randn(B,Number_Rx));
% hh_TR_s=sqrt(1/2)*(randn(B,Number_Rx)+1i*randn(B,Number_Rx));
h_AR_s=sqrt(1/2)*(randn(B,Number_Rx)+1i*randn(B,Number_Rx)); %小尺度衰落
noise0_s= (sqrt(noise_power/2))*((randn(K,B,Number_Rx))+1i*randn(K,B,Number_Rx));
% noise1_s = (sqrt(noise_power/2))*((randn(K,B,Number_Rx))+1i*randn(K,B,Number_Rx));
noise_p_s= (sqrt(noise_power/2))*((randn(K,B,Number_Rx))+1i*randn(K,B,Number_Rx));

theta=-pi+2*pi*rand(K,B,Number_Rx);
tx_signal1_s=cos(theta)+1i*sin(theta);

p_d = zeros(Number_Rx,length(G_vec),B);
p_i_r = zeros(Number_Rx,length(G_vec),B);
p_r = zeros(Number_Rx,length(G_vec),B);
p_i_d = zeros(Number_Rx,length(G_vec),B);
for NoRx=1:Number_Rx
    for i_G=1:length(G_vec)
        GG=G_vec(i_G);
        
        %% CLPC
        theta_CLPC = -pi+2*pi*rand(K,1);
        tx_signal_CLPC = cos(theta_CLPC)+1i*sin(theta_CLPC);
        DTR_CLPC = sqrt((X_Tx-X_Rx_CLPC)^2+(Y_Tx-Y_Rx_CLPC)^2);
        DAR_CLPC = sqrt((X_A-X_Rx_CLPC)^2+(Y_A-Y_Rx_CLPC)^2);
        DTA_CLPC = sqrt((X_A-X_Tx)^2+(Y_A-Y_Tx)^2);
        gpTR_CLPC_dB = -128.1 -37.6*log10(DTR_CLPC/(10^3));
        gpTR_CLPC = (10^(gpTR_CLPC_dB/10));
        gpAR_CLPC_dB = -128.1 -37.6*log10(DAR_CLPC/(10^3));
        gpAR_CLPC = (10^(gpAR_CLPC_dB/10));
        gpTA_CLPC_dB = -128.1 -37.6*log10(DTA_CLPC/(10^3));
        gpTA_CLPC = (10^(gpTA_CLPC_dB/10));
        gsAR_CLPC_dB = normrnd(0,sigma_shadowing);
        gsAR_CLPC = 10^(gsAR_CLPC_dB/10);
        gsTR_CLPC_dB = normrnd(0,sigma_shadowing);
        gsTR_CLPC = 10^(gsTR_CLPC_dB/10);
        gsTA_CLPC_dB = normrnd(0,sigma_shadowing);
        gsTA_CLPC = 10^(gsTA_CLPC_dB/10);
%         h_TR_CLPC = sqrt(1/2)*(randn(1)+1i*randn(1));
        h_TA_CLPC = sqrt(1/2)*(randn(1)+1i*randn(1));
        h_AR_CLPC = sqrt(1/2)*(randn(1)+1i*randn(1));
        noise_0_CLPC = (sqrt(noise_power/2))*((randn(K,1))+1i*randn(K,1));
        noise_p_CLPC = (sqrt(noise_power/2))*((randn(K,1))+1i*randn(K,1));
        
        P0 = gammaT*noise_power/(gpTR_CLPC*gsTR_CLPC);%中继前发射机的发射功率
        
%         %P_D
%         Y_rx_sd = zeros(K,1);
%         for i_k = 1:K
%             Y_rx_sd(i_k)=sqrt(gpTR_CLPC*gsTR_CLPC*P0)*h_TR_CLPC*tx_signal_CLPC(i_k);
%         end
%         P_Y_rx1=Y_rx_sd.*conj(Y_rx_sd);
%         p_d_CLPC = sum(P_Y_rx1)/K;
%         
%         %P_I_R
%         Y_rx_i_r = zeros(K,1);
%         for i_k = 1:K
%             Y_rx_i_r(i_k)=GG*h_TA_CLPC*h_AR_CLPC*sqrt(gpAR_CLPC*gsAR_CLPC*gpTA_CLPC*gsTA_CLPC*P0)*tx_signal_CLPC(i_k)+GG*h_AR_CLPC*sqrt(gpAR_CLPC*gsAR_CLPC)*noise_0_CLPC(i_k)+noise_p_CLPC(i_k); % PR接收的CT的信号和噪声
%         end
%         P_Y_rx_i_r=Y_rx_i_r.*conj(Y_rx_i_r);
%         p_i_r_CLPC=sum(P_Y_rx_i_r)/K;
%         
%         %P_R
%         Y_rx_sr = zeros(K,1);
%         for i_k = 1:K
%             Y_rx_sr(i_k)=GG*h_TA_CLPC*h_AR_CLPC*sqrt(gpTR_CLPC*gsTR_CLPC*gpAR_CLPC*gsAR_CLPC*P0)*tx_signal_CLPC(i_k);   %PR接收的CT的信号（不含噪声）
%         end
%         P_Y_rx2=Y_rx_sr.*conj(Y_rx_sr);   %PR接收的CT的信号（不含噪声）的能量
%         p_r_CLPC=sum(P_Y_rx2)/K;  %PR接收的CT的信号（不含噪声）的功率
% 
%         %P_I_D
%         Y_rx_i_d = zeros(K,1);
%         for i_k = 1:K
%             Y_rx_i_d(i_k)=sqrt(gpTR_CLPC*gsTR_CLPC*P0)*h_TR_CLPC*tx_signal_CLPC(i_k)+GG*h_AR_CLPC*sqrt(gpAR_CLPC*gsAR_CLPC)*noise_0_CLPC(i_k) + noise_p_CLPC(i_k);
%         end
%         P_Y_rx_i_d=Y_rx_i_d.*conj(Y_rx_i_d);
%         p_i_d_CLPC=sum(P_Y_rx_i_d)/K;
        
        PP0 = zeros(K,1);
        for i_k = 1:K
            PP0(i_k)=GG*h_TA_CLPC*h_AR_CLPC*sqrt(gpTA_CLPC*gsTA_CLPC*gpAR_CLPC*gsAR_CLPC*P0)*tx_signal_CLPC(i_k) + sqrt(gpAR_CLPC*gsAR_CLPC*P0)*h_AR_CLPC*tx_signal_CLPC(i_k);  %PR接收的PT和CT的信号的合成（不含噪声）
        end
        PP1=PP0.*conj(PP0);
        PP=sum(PP1)/K;  %PR接收的PT和CT的信号的合成（不含噪声）的功率
        
        P_noi0 = zeros(K,1);
        for i_k = 1:K
            P_noi0(i_k)=GG*h_AR_CLPC*sqrt(gpAR_CLPC*gsAR_CLPC)*noise_0_CLPC(i_k) + noise_p_CLPC(i_k);   %PR接收的PT和CT的噪声
        end
        P_noi1=P_noi0.*conj(P_noi0);
        P_noi=sum(P_noi1)/K;
        gamma_1=PP./P_noi;
        
        % P1
        P1 = gamma_1*noise_power/(gpTR_CLPC*gsTR_CLPC);%CLPC后的发射机功率？？GP GS

        %% Calculate Parameters
        
        for i_B=1:B
            h_TA=h_TA_s(i_B,NoRx);
            h_TR=h_TR_s(i_B,NoRx);
            h_AR=h_AR_s(i_B,NoRx);
            tx_signal0=tx_signal1_s(:,i_B,NoRx);
            noise0= noise0_s(:,i_B,NoRx);
            noise_p= noise_p_s(:,i_B,NoRx);
            
            Y_rx_sd = zeros(K,1);
            for i_k = 1:K
                Y_rx_sd(i_k)=sqrt(gpTR_vec(NoRx)*gsTR_vec(NoRx)*P1)*h_TR*tx_signal0(i_k);
            end
            P_Y_rx1=Y_rx_sd.*conj(Y_rx_sd); 
            p_d(NoRx,i_G,i_B)=sum(P_Y_rx1)/K; %R接收的T直达的信号（不含噪声）的功率

            Y_rx_i_r = zeros(K,1);
            for i_k = 1:K
                Y_rx_i_r(i_k)=GG*h_TA*h_AR*sqrt(gpAR_vec(NoRx)*gsAR_vec(NoRx)*gpTA*gsTA*P1)*tx_signal0(i_k)+GG*h_AR*sqrt(gpAR_vec(NoRx)*gsAR_vec(NoRx))*noise0(i_k)+noise_p(i_k); % PR接收的CT的信号和噪声
            end
            P_Y_rx_i_r=Y_rx_i_r.*conj(Y_rx_i_r); %PR接收的CT信号能量
            p_i_r(NoRx,i_G,i_B)=sum(P_Y_rx_i_r)/K; %(TYD)R接收到T的转发信号和噪声

            Y_rx_sr = zeros(K,1);
            for i_k = 1:K
                Y_rx_sr(i_k)=GG*h_TA*h_AR*sqrt(gpTA*gsTA*gpAR_vec(NoRx)*gsAR_vec(NoRx)*P1)*tx_signal0(i_k);   %PR接收的CT的信号（不含噪声）
            end
            P_Y_rx2=Y_rx_sr.*conj(Y_rx_sr);   %PR接收的CT的信号（不含噪声）的能量
            p_r(NoRx,i_G,i_B)=sum(P_Y_rx2)/K;  %PR接收的CT的信号（不含噪声）的功率

            Y_rx_i_d = zeros(K,1);
            for i_k = 1:K
                Y_rx_i_d(i_k)=sqrt(gpTR_vec(NoRx)*gsTR_vec(NoRx)*P1)*h_TR*tx_signal0(i_k)+GG*h_AR*sqrt(gpAR_vec(NoRx)*gsAR_vec(NoRx))*noise0(i_k)+noise_p(i_k); % PR接收的Tx的信号和噪声
            end
            P_Y_rx_i_d=Y_rx_i_d.*conj(Y_rx_i_d);
            p_i_d(NoRx,i_G,i_B)=sum(P_Y_rx_i_d)/K;
        end
    end
end
p_d_final = zeros(Number_Rx,length(G_vec));
p_i_r_final = zeros(Number_Rx,length(G_vec));
p_r_final = zeros(Number_Rx,length(G_vec));
p_i_d_final = zeros(Number_Rx,length(G_vec));
for NoRx=1:Number_Rx
    for i_G=1:length(G_vec)
        p_d_temp = 0;
        p_i_r_temp = 0;
        p_r_temp = 0;
        p_i_d_temp = 0;
        for i_B=1:B
            p_d_temp = p_d_temp + p_d(NoRx,i_G,i_B);
            p_i_r_temp = p_i_r_temp + p_i_r(NoRx,i_G,i_B);
            p_r_temp = p_r_temp + p_r(NoRx,i_G,i_B);
            p_i_d_temp = p_i_d_temp + p_i_d(NoRx,i_G,i_B);
        end
        p_d_final(NoRx,i_G) = p_d_temp/B;
        p_i_r_final(NoRx,i_G) = p_i_r_temp/B;
        p_r_final(NoRx,i_G) = p_r_temp/B;
        p_i_d_final(NoRx,i_G) = p_i_d_temp/B;
    end
end
if(Time_Measure)
    t = toc;
    fprintf('Channel Analysis Time = %i\n',t);
end
end

