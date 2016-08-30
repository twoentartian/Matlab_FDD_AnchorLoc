function [ X_FinalPoint,Y_FinalPoint,X_Rx_Success,Y_Rx_Success ] = FindByNetMethod( MapLength,Number_Rx,p_i_d_final,p_d_final,p_r_final,p_i_r_final,Times_From_A,Times_From_Tx,Threshold_Time,X_Tx,Y_Tx,X_Rx,Y_Rx,i_G )
%FINDBYNETMETHOD (MapLength,Number_Rx,p_i_d_final,p_d_final,p_r_final,p_i_r_final,Times_From_A,Times_From_Tx,Threshold_Time,X_Tx,Y_Tx,X_Rx,Y_Rx,i_G)
%   此处显示详细说明

% Const
c = 3*10^8;

% Check arguements
AllSteps = 50;
Power_Addition = 2;% it must be an EVEN number

% Precise Iteration arguements
Precision_Iteration_Time = 20;
Precision_Iteration_AllSteps = 10;

%% Success Determine
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

%% Rough Iteration
X_MaybePoints = zeros(AllSteps);
Y_MaybePoints = zeros(AllSteps);
% Generate the Location of Find Points
for x = 1:AllSteps
    for y = 1:AllSteps
        X_MaybePoints(x,y) = MapLength/(AllSteps-1)*(x-1)*2 - MapLength;
        Y_MaybePoints(x,y) = MapLength/(AllSteps-1)*(y-1)*2 - MapLength;
    end
end
% Find the Minimum Value
Min_Sign_Value = Inf;
X_Maybe_Min = Inf;
Y_Maybe_Min = Inf;
Sign_Value = zeros(AllSteps,AllSteps);
for x = 1:AllSteps
    for y = 1:AllSteps
        for NoRx = 1:Number_Rx
            if p_d_final(NoRx,i_G)/p_i_r_final(NoRx,i_G) < p_r_final(NoRx,i_G)/p_i_d_final(NoRx,i_G) && Threshold_Time < abs(Times_From_A(NoRx) - Times_From_Tx(NoRx))
                Time_Temp = (sqrt((X_Tx-X_MaybePoints(x,y))^2 + (Y_Tx-Y_MaybePoints(x,y))^2) + sqrt((X_Rx(NoRx)-X_MaybePoints(x,y))^2 + (Y_Rx(NoRx)-Y_MaybePoints(x,y))^2))/c;
                Sign_Value(x,y) = Sign_Value(x,y) + (Time_Temp*10^7 - Times_From_A(NoRx)*10^7)^Power_Addition;
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
if X_Maybe_Min==Inf && Y_Maybe_Min==Inf
    fprintf('Success Ratio is Zero, Failed to locate.\n');
    X_FinalPoint = Inf;
    Y_FinalPoint = Inf;
    return;
end
Last_Iteration_Step = MapLength/(AllSteps-1)*2;%the step of base finding
X_FinalPoint = X_MaybePoints(X_Maybe_Min,Y_Maybe_Min);
Y_FinalPoint = Y_MaybePoints(X_Maybe_Min,Y_Maybe_Min);

%% Precise Iteration (Optional)
for Iteration_Time = 1:Precision_Iteration_Time
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
                    Sign_Value_Iteration(x,y) = Sign_Value_Iteration(x,y) + (Time_Temp*10^7 - Times_From_A(NoRx)*10^7)^Power_Addition;
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

end

