function [ X_FinalPoint,Y_FinalPoint,X_Rx_Success,Y_Rx_Success ] = FindByDiffMethod( MapLength,Number_Rx,p_i_d_final,p_d_final,p_r_final,p_i_r_final,Times_From_A,Times_From_Tx,Threshold_Time,X_Tx,Y_Tx,X_Rx,Y_Rx,i_G )
%FINDBYDIFFMETHOD 此处显示有关此函数的摘要
%   此处显示详细说明

% Diff:dF/dX
% 2*((2*X_A - 2*X_Tx)/(2*((X_A - X_Tx)^2 + (Y_A - Y_Tx)^2)^(1/2)) + (2*X_A - 2*X_Rx_Success)/(2*((X_A - X_Rx_Success)^2 + (Y_A - Y_Rx_Success)^2)^(1/2)))*(((X_A - X_Tx)^2 + (Y_A - Y_Tx)^2)^(1/2) - Time*c + ((X_A - X_Rx_Success)^2 + (Y_A - Y_Rx_Success)^2)^(1/2))
% Diff:dF/dY
% 2*((2*Y_A - 2*Y_Tx)/(2*((X_A - X_Tx)^2 + (Y_A - Y_Tx)^2)^(1/2)) + (2*Y_A - 2*Y_Rx_Success)/(2*((X_A - X_Rx_Success)^2 + (Y_A - Y_Rx_Success)^2)^(1/2)))*(((X_A - X_Tx)^2 + (Y_A - Y_Tx)^2)^(1/2) - Time*c + ((X_A - X_Rx_Success)^2 + (Y_A - Y_Rx_Success)^2)^(1/2))

% Const
c = 3*10^8;

% Check arguements
Iteration_Time = 1;
Iteration_AllSteps = 100;

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


%% Main
X_FinalPoint = 0;
Y_FinalPoint = 0;
Last_Iteration_Step = MapLength;%the step of base finding
for i_Iteration_Time = 1:Iteration_Time
    X_MaybePoints = zeros(Iteration_AllSteps);
    Y_MaybePoints = zeros(Iteration_AllSteps);
    for x = 1:Iteration_AllSteps
        for y = 1:Iteration_AllSteps
            X_MaybePoints(x,y) = Last_Iteration_Step/(Iteration_AllSteps-1)*(x-1)*2 - Last_Iteration_Step + X_FinalPoint;
            Y_MaybePoints(x,y) = Last_Iteration_Step/(Iteration_AllSteps-1)*(y-1)*2 - Last_Iteration_Step + Y_FinalPoint;
        end
    end
    
    Sign_Min = Inf;
    SignValue = zeros(Iteration_AllSteps);
    X_SignValue = zeros(Iteration_AllSteps);
    Y_SignValue = zeros(Iteration_AllSteps);
    
    for x = 1:Iteration_AllSteps
        for y = 1:Iteration_AllSteps
            for NoRx = 1:SuccessCounter
                Time = Times_From_A(NoRx);
                % f = @(X_A,Y_A) (((X_A-X_Rx_Success(NoRx))^2+(Y_A-Y_Rx_Success(NoRx))^2)^0.5 + ((X_A-X_Tx)^2+(Y_A-Y_Tx)^2)^0.5 - Time * c)^2;
                % dfx = @(X_A,Y_A) 2*((2*X_A - 2*X_Tx)/(2*((X_A - X_Tx)^2 + (Y_A - Y_Tx)^2)^(1/2)) + (2*X_A - 2*X_Rx_Success(NoRx))/(2*((X_A - X_Rx_Success(NoRx))^2 + (Y_A - Y_Rx_Success(NoRx))^2)^(1/2)))*(((X_A - X_Tx)^2 + (Y_A - Y_Tx)^2)^(1/2) - Time*c + ((X_A - X_Rx_Success(NoRx))^2 + (Y_A - Y_Rx_Success(NoRx))^2)^(1/2));
                % dfy = @(X_A,Y_A) 2*((2*Y_A - 2*Y_Tx)/(2*((X_A - X_Tx)^2 + (Y_A - Y_Tx)^2)^(1/2)) + (2*Y_A - 2*Y_Rx_Success(NoRx))/(2*((X_A - X_Rx_Success(NoRx))^2 + (Y_A - Y_Rx_Success(NoRx))^2)^(1/2)))*(((X_A - X_Tx)^2 + (Y_A - Y_Tx)^2)^(1/2) - Time*c + ((X_A - X_Rx_Success(NoRx))^2 + (Y_A - Y_Rx_Success(NoRx))^2)^(1/2));
                temp_Diff_X = dfx(X_MaybePoints(x,y),Y_MaybePoints(x,y),X_Tx,Y_Tx,X_Rx_Success,Y_Rx_Success,NoRx,Time,c);
                temp_Diff_Y = dfy(X_MaybePoints(x,y),Y_MaybePoints(x,y),X_Tx,Y_Tx,X_Rx_Success,Y_Rx_Success,NoRx,Time,c);
                X_SignValue(x,y) = X_SignValue(x,y) + temp_Diff_X;
                Y_SignValue(x,y) = Y_SignValue(x,y) + temp_Diff_Y;
                SignValue(x,y) = X_SignValue(x,y) + Y_SignValue(x,y);
            end

            if(abs(SignValue(x,y)) < Sign_Min)
                X_Min = X_MaybePoints(x,y);
                Y_Min = Y_MaybePoints(x,y);
                Sign_Min = abs(SignValue(x,y));
            end
        end
    end
    
    %figure();
    %mesh(X_SignValue);
    %figure();
    %mesh(Y_SignValue);
    
    X_FinalPoint = X_Min;
    Y_FinalPoint = Y_Min;
    Last_Iteration_Step = Last_Iteration_Step*2/(Iteration_AllSteps-1)*2;%the step of Iteration
    
    fprintf('X=%f Y=%f\n',X_FinalPoint,Y_FinalPoint);
end
disp(X_Check(X_Tx,Y_Tx,X_Rx_Success,Y_Rx_Success,Time,c,SuccessCounter));
X = linspace(-100,100);
Y = linspace(-100,100);
Z = zeros(length(X),length(Y));
for i_X = 1:length(X)
    for i_Y = 1:length(Y)
        for NoRx = 1:SuccessCounter
            Z(i_X,i_Y) = Z(i_X,i_Y) + f(X(i_X),Y(i_Y),X_Tx,Y_Tx,X_Rx_Success,Y_Rx_Success,NoRx,Time,c);
        end
    end
end
figure();
mesh(X,Y,Z);
xlabel('X');
ylabel('Y');

end

function Out = f(X_A,Y_A,X_Tx,Y_Tx,X_Rx_Success,Y_Rx_Success,NoRx,Time,c)
Out = abs(sqrt((X_A-X_Rx_Success(NoRx))^2+(Y_A-Y_Rx_Success(NoRx))^2) + sqrt((X_A-X_Tx)^2+(Y_A-Y_Tx)^2) - Time * c);
end

function Out = Check_f(X_A,Y_A,X_Tx,Y_Tx,X_Rx_Success,Y_Rx_Success,Time,c,SuccessCounter)
Out = 0;
for NoRx = 1:SuccessCounter
    Out = Out + f(X_A,Y_A,X_Tx,Y_Tx,X_Rx_Success,Y_Rx_Success,NoRx,Time,c);
end
end

function Out = dfx(X_A,Y_A,X_Tx,Y_Tx,X_Rx_Success,Y_Rx_Success,NoRx,Time,c)
Out = 2*((2*X_A - 2*X_Tx)/(2*((X_A - X_Tx)^2 + (Y_A - Y_Tx)^2)) + (2*X_A - 2*X_Rx_Success(NoRx))/(2*((X_A - X_Rx_Success(NoRx))^2 + (Y_A - Y_Rx_Success(NoRx))^2)^(1/2)))*(((X_A - X_Tx)^2 + (Y_A - Y_Tx)^2)^(1/2) - Time*c + ((X_A - X_Rx_Success(NoRx))^2 + (Y_A - Y_Rx_Success(NoRx))^2)^(1/2));
end

function Out = dfy(X_A,Y_A,X_Tx,Y_Tx,X_Rx_Success,Y_Rx_Success,NoRx,Time,c)
Out = 2*((2*Y_A - 2*Y_Tx)/(2*((X_A - X_Tx)^2 + (Y_A - Y_Tx)^2)^(1/2)) + (2*Y_A - 2*Y_Rx_Success(NoRx))/(2*((X_A - X_Rx_Success(NoRx))^2 + (Y_A - Y_Rx_Success(NoRx))^2)^(1/2)))*(((X_A - X_Tx)^2 + (Y_A - Y_Tx)^2)^(1/2) - Time*c + ((X_A - X_Rx_Success(NoRx))^2 + (Y_A - Y_Rx_Success(NoRx))^2)^(1/2));
end

function Out = X_Check(X_Tx,Y_Tx,X_Rx_Success,Y_Rx_Success,Time,c,SuccessCounter)
Out = 0;
for NoRx = 1:SuccessCounter
    Out = Out + 2*((2*80 - 2*X_Tx)/(2*((80 - X_Tx)^2 + (80 - Y_Tx)^2)^(1/2)) + (2*80 - 2*X_Rx_Success(NoRx))/(2*((80 - X_Rx_Success(NoRx))^2 + (80 - Y_Rx_Success(NoRx))^2)^(1/2)))*(((80 - X_Tx)^2 + (80 - Y_Tx)^2)^(1/2) - Time*c + ((80 - X_Rx_Success(NoRx))^2 + (80 - Y_Rx_Success(NoRx))^2)^(1/2));
end
end

function Out = X_Diff(X_A,Y_A,X_Tx,Y_Tx,X_Rx_Success,Y_Rx_Success,NoRx,Time,c)
Delta = 0.000000001;
f = @(X_A,Y_A) (((X_A-X_Rx_Success(NoRx))^2+(Y_A-Y_Rx_Success(NoRx))^2)^0.5 + ((X_A-X_Tx)^2+(Y_A-Y_Tx)^2)^0.5 - Time * c)^2;
Out = (f(X_A+Delta,Y_A) - f(X_A,Y_A))/Delta;
end