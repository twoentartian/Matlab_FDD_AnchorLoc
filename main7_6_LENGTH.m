%%
clearvars;
clc;

%% Constant
c = 3*10^8;

%% Init
% Map
MapLength_vec = [100 500 1000];
sigma_shadowing=4;
B=300;
K=300;
gammaT_dB=10;
gammaT=10^(gammaT_dB/10);
noise_power_dB=-144;
noise_power=10^(noise_power_dB/10);

Number_Map = 10;

% Check arguements
AllSteps = 50;
Power_Addition = 2;% it must be an EVEN number

% Precise Iteration arguements
Precision_Iteration_Time = 20;
Precision_Iteration_AllSteps = 10;

% Anchor
X_A = -20;
Y_A = 0;
Delay_Time_Dis = 5;
Delay_Time = Delay_Time_Dis/c;
Delay_Time_k = 0.1;

G_dB_vec = 15:5:40;%Power Gain
G_vec = 10.^(G_dB_vec/10);

% Rx(CLPC)
X_Rx_CLPC = 0;
Y_Rx_CLPC = 0;

% Rx(slience)
Number_Rx = 100;
X_Rx = zeros(Number_Rx,1);
Y_Rx = zeros(Number_Rx,1);
Times_From_A = zeros(Number_Rx,1);
Times_From_Tx = zeros(Number_Rx,1);

% TX
X_Tx = 20;
Y_Tx = 0;

% Success Counter
CounterSuccess = 0;

Threshold_Distance = 1.5;
Threshold_Time = Threshold_Distance / c; % 5ns

SuccessRatio_final_vec = zeros(length(MapLength_vec),length(G_dB_vec));
RMSE_final_vec = Inf(length(MapLength_vec),length(G_dB_vec));
for i_Length = 1:1:length(MapLength_vec)
    MapLength = MapLength_vec(i_Length);
    
	%% MAIN
	RMSE_vec = Inf(Number_Map,length(G_dB_vec));
	SuccessRatio_vec = zeros(Number_Map,length(G_dB_vec));
	for i_Map = 1:Number_Map
		%% Generate Map
		for NoRx = 1:Number_Rx
			X_Rx(NoRx) = rand()*2*MapLength - MapLength;
			Y_Rx(NoRx) = rand()*2*MapLength - MapLength;
			Times_From_Tx(NoRx) = sqrt((X_Tx-X_Rx(NoRx))^2 + (Y_Tx-Y_Rx(NoRx))^2)/c;
			Times_From_A(NoRx) = (sqrt((X_Tx-X_A)^2 + (Y_Tx-Y_A)^2) + sqrt((X_A-X_Rx(NoRx))^2 + (Y_A-Y_Rx(NoRx))^2))/c + Delay_Time*(1+randn(1)*Delay_Time_k);
		end
		
		%% Generation Parameters
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
		hh_TA_s=sqrt(1/2)*(randn(B,Number_Rx)+1i*randn(B,Number_Rx));
		h_TR_s=sqrt(1/2)*(randn(B,Number_Rx)+1i*randn(B,Number_Rx));
		hh_TR_s=sqrt(1/2)*(randn(B,Number_Rx)+1i*randn(B,Number_Rx));
		h_AR_s=sqrt(1/2)*(randn(B,Number_Rx)+1i*randn(B,Number_Rx)); %小尺度衰落
		noise0_s= (sqrt(noise_power/2))*((randn(K,B,Number_Rx))+1i*randn(K,B,Number_Rx));
		noise1_s = (sqrt(noise_power/2))*((randn(K,B,Number_Rx))+1i*randn(K,B,Number_Rx));
		noise_p_s= (sqrt(noise_power/2))*((randn(K,B,Number_Rx))+1i*randn(K,B,Number_Rx));
		
		theta1=-pi+2*pi*rand(K,B,Number_Rx);
		tx_signal0_s=cos(theta1)+1i*sin(theta1);
		
		p_d = zeros(Number_Rx,length(G_vec),B);
		p_i_r = zeros(Number_Rx,length(G_vec),B);
		p_r = zeros(Number_Rx,length(G_vec),B);
		p_i_d = zeros(Number_Rx,length(G_vec),B);
		for NoRx=1:Number_Rx
			for i_G=1:length(G_vec)
				GG=G_vec(i_G);
                
                % CLPC
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
                h_TR_CLPC = sqrt(1/2)*(randn(1)+1i*randn(1));
                h_TA_CLPC = sqrt(1/2)*(randn(1)+1i*randn(1));
                h_AR_CLPC = sqrt(1/2)*(randn(1)+1i*randn(1));
                noise_0_CLPC = (sqrt(noise_power/2))*((randn(K,1))+1i*randn(K,1));
                noise_p_CLPC = (sqrt(noise_power/2))*((randn(K,1))+1i*randn(K,1));
            
                P0 = gammaT*noise_power/(gpTR_CLPC*gsTR_CLPC);%中继前发射机的发射功率
            
                %P_D
                Y_rx_sd = zeros(K,1);
                for i_k = 1:K
                    Y_rx_sd(i_k)=sqrt(gpTR_CLPC*gsTR_CLPC*P0)*h_TR_CLPC*tx_signal_CLPC(i_k);
                end
                P_Y_rx1=Y_rx_sd.*conj(Y_rx_sd);
                p_d_CLPC = sum(P_Y_rx1)/K;
            
                %P_I_R
                Y_rx_i_r = zeros(K,1);
                for i_k = 1:K
                    Y_rx_i_r(i_k)=GG*h_TA_CLPC*h_AR_CLPC*sqrt(gpAR_CLPC*gsAR_CLPC*gpTA_CLPC*gsTA_CLPC*P0)*tx_signal_CLPC(i_k)+GG*h_AR_CLPC*sqrt(gpAR_CLPC*gsAR_CLPC)*noise_0_CLPC(i_k)+noise_p_CLPC(i_k); % PR接收的CT的信号和噪声
                end
                P_Y_rx_i_r=Y_rx_i_r.*conj(Y_rx_i_r);
                p_i_r_CLPC=sum(P_Y_rx_i_r)/K;
            
                %P_R
                Y_rx_sr = zeros(K,1);
                for i_k = 1:K
                    Y_rx_sr(i_k)=GG*h_TA_CLPC*h_AR_CLPC*sqrt(gpTR_CLPC*gsTR_CLPC*gpAR_CLPC*gsAR_CLPC*P0)*tx_signal_CLPC(i_k);   %PR接收的CT的信号（不含噪声）
                end
                P_Y_rx2=Y_rx_sr.*conj(Y_rx_sr);   %PR接收的CT的信号（不含噪声）的能量
                p_r_CLPC=sum(P_Y_rx2)/K;  %PR接收的CT的信号（不含噪声）的功率
            
                %P_I_D
                Y_rx_i_d = zeros(K,1);
                for i_k = 1:K
                    Y_rx_i_d(i_k)=sqrt(gpTR_CLPC*gsTR_CLPC*P0)*h_TR_CLPC*tx_signal_CLPC(i_k)+GG*h_AR_CLPC*sqrt(gpAR_CLPC*gsAR_CLPC)*noise_0_CLPC(i_k) + noise_p_CLPC(i_k);
                end
                P_Y_rx_i_d=Y_rx_i_d.*conj(Y_rx_i_d);
                p_i_d_CLPC=sum(P_Y_rx_i_d)/K;
                
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
	
				for i_B=1:B
					h_TA=h_TA_s(i_B,NoRx);
					h_TR=h_TR_s(i_B,NoRx);
					h_AR=h_AR_s(i_B,NoRx);
					tx_signal0=tx_signal0_s(:,i_B,NoRx);
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
		
		for i_G = 1:length(G_dB_vec)
			GG = G_dB_vec(i_G);
			X_Rx_Success = Inf(Number_Rx,1);
			Y_Rx_Success = Inf(Number_Rx,1);
			SuccessCounter = 0;
            for NoRx=1:Number_Rx
				if p_d_final(NoRx,i_G)/p_i_r_final(NoRx,i_G) < p_r_final(NoRx,i_G)/p_i_d_final(NoRx,i_G) && Threshold_Time < abs(Times_From_A(NoRx) - Times_From_Tx(NoRx))
					X_Rx_Success(SuccessCounter + 1) = X_Rx(NoRx);
					Y_Rx_Success(SuccessCounter + 1) = Y_Rx(NoRx);
					SuccessCounter = SuccessCounter + 1;
				end
            end
			
            %% Main Process
			fprintf('-----LOG-----\n');
			fprintf('Map: %i   Anchor Power: %idB\n',i_Map,GG);
			
			X_MaybePoints = zeros(AllSteps);
			Y_MaybePoints = zeros(AllSteps);
			% Generate the Location of Find Points
			fprintf('Generate the Location of Net.\n');
			for x = 1:AllSteps
				for y = 1:AllSteps
					X_MaybePoints(x,y) = MapLength/(AllSteps-1)*(x-1)*2 - MapLength;
					Y_MaybePoints(x,y) = MapLength/(AllSteps-1)*(y-1)*2 - MapLength;
				end
			end
			% Find the Minimum Value
			fprintf('Find the Minimum Value\n');
			Min_Sign_Value = Inf;
			X_Maybe_Min = Inf;
			Y_Maybe_Min = Inf;
			Sign_Value = zeros(AllSteps,AllSteps);
			for x = 1:AllSteps
				for y = 1:AllSteps
					for NoRx = 1:Number_Rx
						if p_d_final(NoRx,i_G)/p_i_r_final(NoRx,i_G) < p_r_final(NoRx,i_G)/p_i_d_final(NoRx,i_G) && Threshold_Time < abs(Times_From_A(NoRx) - Times_From_Tx(NoRx))
							Time_Temp = (sqrt((X_Tx-X_MaybePoints(x,y))^2 + (Y_Tx-Y_MaybePoints(x,y))^2) + sqrt((X_Rx(NoRx)-X_MaybePoints(x,y))^2 + (Y_Rx(NoRx)-Y_MaybePoints(x,y))^2))/c;
							Sign_Value(x,y) = Sign_Value(x,y) + (Time_Temp*10^7 - (Times_From_A(NoRx)-Delay_Time)*10^7)^Power_Addition;
						end
					end
					Sign_Value(x,y) = Sign_Value(x,y)/SuccessCounter;
					if Sign_Value(x,y) < Min_Sign_Value
						Min_Sign_Value = Sign_Value(x,y);
						X_Maybe_Min = x;
						Y_Maybe_Min = y;
					end
				end
			end
			if X_Maybe_Min==Inf && Y_Maybe_Min==Inf && CounterSuccess==0
				fprintf('Success Ratio is Zero, Failed to locate.\n');
				continue;
			end
			Last_Iteration_Step = MapLength/(AllSteps-1)*2;%the step of base finding
			X_FinalPoint = X_MaybePoints(X_Maybe_Min,Y_Maybe_Min);
			Y_FinalPoint = Y_MaybePoints(X_Maybe_Min,Y_Maybe_Min);
			
           %% Precise Iteration (Optional)
			for Iteration_Time = 1:Precision_Iteration_Time
				%fprintf('Precise Iteration for %i\n',Iteration_Time);
				
				% Generate the Location of Find Points
				X_MaybePoints_Iteration = zeros(Precision_Iteration_AllSteps);
				Y_MaybePoints_Iteration = zeros(Precision_Iteration_AllSteps);
				for x = 1:Precision_Iteration_AllSteps
					for y = 1:Precision_Iteration_AllSteps
						X_MaybePoints_Iteration(x,y) = Last_Iteration_Step/(Precision_Iteration_AllSteps-1)*(x-1)*2 - Last_Iteration_Step + X_FinalPoint;
						Y_MaybePoints_Iteration(x,y) = Last_Iteration_Step/(Precision_Iteration_AllSteps-1)*(y-1)*2 - Last_Iteration_Step + Y_FinalPoint;
					end
				end
				
				% Find the Minimum Value
				Min_Sign_Value = Inf;
				X_Maybe_Min = Inf;
				Y_Maybe_Min = Inf;
				Sign_Value_Iteration = zeros(Precision_Iteration_AllSteps,Precision_Iteration_AllSteps);
				for x = 1:Precision_Iteration_AllSteps
					for y = 1:Precision_Iteration_AllSteps
						for NoRx = 1:Number_Rx
							if p_d_final(NoRx,i_G)/p_i_r_final(NoRx,i_G) < p_r_final(NoRx,i_G)/p_i_d_final(NoRx,i_G) && Threshold_Time < abs(Times_From_A(NoRx) - Times_From_Tx(NoRx))
								Time_Temp = (sqrt((X_Tx-X_MaybePoints_Iteration(x,y))^2 + (Y_Tx-Y_MaybePoints_Iteration(x,y))^2) + sqrt((X_Rx(NoRx)-X_MaybePoints_Iteration(x,y))^2 + (Y_Rx(NoRx)-Y_MaybePoints_Iteration(x,y))^2))/c;
								Sign_Value_Iteration(x,y) = Sign_Value_Iteration(x,y) + (Time_Temp*10^7 - (Times_From_A(NoRx)-Delay_Time)*10^7)^Power_Addition;
							end
						end
						Sign_Value_Iteration(x,y) = Sign_Value_Iteration(x,y)/SuccessCounter;
						if Sign_Value_Iteration(x,y) < Min_Sign_Value
							Min_Sign_Value = Sign_Value_Iteration(x,y);
							X_Maybe_Min = x;
							Y_Maybe_Min = y;
						end
					end
				end
				if X_Maybe_Min==Inf && Y_Maybe_Min==Inf
					error('Impossable math situation!');
				end
				Last_Iteration_Step = Last_Iteration_Step*2/(Precision_Iteration_AllSteps-1)*2;%the step of Iteration
				X_FinalPoint = X_MaybePoints_Iteration(X_Maybe_Min,Y_Maybe_Min);
				Y_FinalPoint = Y_MaybePoints_Iteration(X_Maybe_Min,Y_Maybe_Min);
			end
			RMSE_vec(i_Map,i_G) = (X_FinalPoint - X_A)^2 + (Y_FinalPoint - Y_A)^2;
			if(SuccessCounter >= 3)
				SuccessRatio_vec(i_Map,i_G) = 1;
			else
				SuccessRatio_vec(i_Map,i_G) = 0;
			end
		end
	end
	%% Plot
    for i_G = 1:length(G_dB_vec)
		temp = sum(SuccessRatio_vec,1);
		SuccessRatio_final_vec(i_Length,i_G) = temp(i_G)/Number_Map;
		
		temp = sum(RMSE_vec,1);
		if(SuccessRatio_final_vec(i_Length,i_G) == 0)
			RMSE_final_vec(i_Length,i_G) = Inf;
		else
			counter = 0;
			temp = 0;
			for i_Map = 1:Number_Map
				if(RMSE_vec(i_Map,i_G) ~= Inf)
					counter = counter+1;
					temp = temp + RMSE_vec(i_Map,i_G);
				else
					
				end
			end
			RMSE_final_vec(i_Length,i_G) = sqrt(temp/counter);
		end
    end
end

figure();
hold on;
grid on;
for i_Length = 1:1:length(MapLength_vec)
    plot(G_dB_vec,RMSE_final_vec(i_Length,:),'-');
end
xlabel('Gain(dB)');
axis([min(G_dB_vec) max(G_dB_vec) -Inf Inf]);
ylabel('RMSE');

figure();
hold on;
grid on;
for i_Length = 1:1:length(MapLength_vec)
    plot(G_dB_vec,SuccessRatio_final_vec(i_Length,:),'-');
end
axis([min(G_dB_vec) max(G_dB_vec) -Inf Inf]);
xlabel('Gain(dB)');
ylabel('Success Ratio');