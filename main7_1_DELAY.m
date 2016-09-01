%%
clearvars;
clc;

%% Constant
c = 3*10^8;

%% Init
% Map
MapLength = 100;
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

G_dB_vec = 25:1:35;%Power Gain
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

% Rx Parameters
Threshold_Distance = 1.5;
Threshold_Time = Threshold_Distance / c; % 5ns

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
    [ p_d_final,p_i_r_final,p_r_final,p_i_d_final ] = Channel_Analusis( X_Tx,Y_Tx,X_Rx,Y_Rx,X_A,Y_A,Number_Rx,MapLength,G_vec );
	
    for i_G = 1:length(G_dB_vec)
        fprintf('Gain = %i Map = %i\n',G_dB_vec(i_G),i_Map);
		GG = G_dB_vec(i_G);
        [X_FinalPoint,Y_FinalPoint,Success_Set,SuccessCounter] = FindByNetMethod( MapLength,Number_Rx,p_i_d_final,p_d_final,p_r_final,p_i_r_final,Times_From_A,Times_From_Tx,Threshold_Time,X_Tx,Y_Tx,X_Rx,Y_Rx,i_G );
        %[X_FinalPoint,Y_FinalPoint,Success_Set,SuccessCounter] = FindByDiffMethod( MapLength,Number_Rx,p_i_d_final,p_d_final,p_r_final,p_i_r_final,Times_From_A,Times_From_Tx,Threshold_Time,X_Tx,Y_Tx,X_Rx,Y_Rx,i_G );

        RMSE_vec(i_Map,i_G) = sqrt((X_FinalPoint - X_A)^2 + (Y_FinalPoint - Y_A)^2);
        if(SuccessCounter >= 3)
            SuccessRatio_vec(i_Map,i_G) = 1;
        else
            SuccessRatio_vec(i_Map,i_G) = 0;
        end
    end
end
%% Plot
RMSE_final_vec = Inf(length(G_dB_vec),1);
SuccessRatio_final_vec = zeros(length(G_dB_vec),1);
for i_G = 1:length(G_dB_vec)
    temp = sum(SuccessRatio_vec,1);
    SuccessRatio_final_vec(i_G) = temp(i_G)/Number_Map;
    
    temp = sum(RMSE_vec,1);
    if(SuccessRatio_final_vec(i_G) == 0)
        RMSE_final_vec(i_G) = Inf;
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
        RMSE_final_vec(i_G) = temp/counter;
    end
end
figure();
hold on;
grid on;
plot(G_dB_vec,RMSE_final_vec,'-');
axis([min(G_dB_vec) max(G_dB_vec) -Inf Inf]);
xlabel('Gain(dB)');
ylabel('RMSE');

figure();
hold on;
grid on;
plot(G_dB_vec,SuccessRatio_final_vec,'-');
axis([min(G_dB_vec) max(G_dB_vec) -Inf Inf]);
xlabel('Gain(dB)');
ylabel('Success Ratio');