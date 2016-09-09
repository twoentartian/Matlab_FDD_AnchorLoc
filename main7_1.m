%%
clearvars;
clc;

%% Constant
c = 3*10^8;

%% Init
% Map
MapLength = 100;

% Anchor
X_A = rand()*2*MapLength - MapLength;
Y_A = rand()*2*MapLength - MapLength;

%X_A = 99.1;
%Y_A = 99.1;

G_dB_vec = 80:5:81;%Power Gain
G_vec = 10.^(G_dB_vec/10);

% Rx(slience)
Number_Rx = 100;
X_Rx = zeros(Number_Rx,1);
Y_Rx = zeros(Number_Rx,1);
Times_From_A = zeros(Number_Rx,1);
Times_From_Tx = zeros(Number_Rx,1);

% Min Time(RX)
Threshold_Distance = 1.5;
Threshold_Time = Threshold_Distance / c; % 5ns

% TX
X_Tx = 0;
Y_Tx = 0;

%% Generate Map
for NoRx = 1:Number_Rx
    X_Rx(NoRx) = rand()*2*MapLength - MapLength;
    Y_Rx(NoRx) = rand()*2*MapLength - MapLength;
    Times_From_Tx(NoRx) = sqrt((X_Tx-X_Rx(NoRx))^2 + (Y_Tx-Y_Rx(NoRx))^2)/c;
    Times_From_A(NoRx) = (sqrt((X_Tx-X_A)^2 + (Y_Tx-Y_A)^2) + sqrt((X_A-X_Rx(NoRx))^2 + (Y_A-Y_Rx(NoRx))^2))/c;
end

%% Generat Parameters
[ p_d_final,p_i_r_final,p_r_final,p_i_d_final ] = Channel_Analusis( X_Tx,Y_Tx,X_Rx,Y_Rx,X_A,Y_A,Number_Rx,MapLength,G_vec );

for i_G = 1:length(G_dB_vec)
	GG = G_dB_vec(i_G);
    
    fprintf('Anchor Power: %ddB\n',GG);
    [ X_FinalPoint,Y_FinalPoint,Success_Set,SuccessCounter,Time ] = FindByNetMethod( MapLength,Number_Rx,p_i_d_final,p_d_final,p_r_final,p_i_r_final,Times_From_A,Times_From_Tx,Threshold_Time,X_Tx,Y_Tx,X_Rx,Y_Rx,i_G );
    [ X_FinalPoint,Y_FinalPoint,Success_Set,SuccessCounter,Time ] = FindByDiffMethod( MapLength,Number_Rx,p_i_d_final,p_d_final,p_r_final,p_i_r_final,Times_From_A,Times_From_Tx,Threshold_Time,X_Tx,Y_Tx,X_Rx,Y_Rx,i_G );
    
    %% Plot
	figure();
	hold on;
	axis([-MapLength MapLength -MapLength MapLength]);
    title(sprintf('Gain=%ddB Final Point:(%d,%d)',GG,round(X_FinalPoint),round(Y_FinalPoint)));
	plot(X_Tx,Y_Tx,'ob');
	plot(X_A,Y_A,'og');
	plot(X_Rx,Y_Rx,'oc');
    for i = 1:SuccessCounter
        plot(X_Rx(Success_Set(i)),Y_Rx(Success_Set(i)),'ok');
    end
	plot(X_FinalPoint,Y_FinalPoint,'or');
	grid on;
end

%% Save
%save data; 
	