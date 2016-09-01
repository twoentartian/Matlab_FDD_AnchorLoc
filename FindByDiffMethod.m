function [ X_FinalPoint,Y_FinalPoint,Success_Set,SuccessCounter ] = FindByDiffMethod( MapLength,Number_Rx,p_i_d_final,p_d_final,p_r_final,p_i_r_final,Times_From_A,Times_From_Tx,Threshold_Time,X_Tx,Y_Tx,X_Rx,Y_Rx,i_G )
%FINDBYDIFFMETHOD Typical time use:xe-02 s or xe-03 s
%   此处显示详细说明
Time_Measure = true;
if(Time_Measure)
    tic;
end

% Const
c = 3*10^8;
Error = 0.0000000000001;
IterationTime = 1;

%% Success Determine
Success_Set = Inf(Number_Rx,1);
SuccessCounter = 0;
for NoRx=1:Number_Rx
    if p_d_final(NoRx,i_G)/p_i_r_final(NoRx,i_G) < p_r_final(NoRx,i_G)/p_i_d_final(NoRx,i_G) && Threshold_Time < abs(Times_From_A(NoRx) - Times_From_Tx(NoRx))
        SuccessCounter = SuccessCounter + 1;
        Success_Set(SuccessCounter) = NoRx;
    end
end

if SuccessCounter==0
    fprintf('Success Ratio is Zero, Failed to locate.\n');
    X_FinalPoint = Inf;
    Y_FinalPoint = Inf;
    return;
end

%% Main
% x = -100:1:100;
% y = -100:1:100;
% z = -Inf(length(x),length(y));
% for xTemp = 1:length(x)
%     for yTemp = 1:length(y)
%         z(xTemp,yTemp) = f_Success(x(xTemp),y(yTemp),X_Tx,Y_Tx,X_Rx,Y_Rx,Success_Set,SuccessCounter,Times_From_A,c);
%     end
% end
% 
% mesh(z);

X_FinalResult = zeros(IterationTime,1);
Y_FinalResult = zeros(IterationTime,1);
for time = 1:IterationTime
    while(true)
        X_FinalPoint = rand()*MapLength*2 - MapLength;
        Y_FinalPoint = rand()*MapLength*2 - MapLength;
        SuccessSign = false;
        while(true)
            Old_X_FinalPoint = X_FinalPoint;
            Old_Y_FinalPoint = Y_FinalPoint;
            [dfx,dfy] = find_dfx_dfy(X_FinalPoint,Y_FinalPoint,X_Tx,Y_Tx,X_Rx,Y_Rx,Success_Set,SuccessCounter,Times_From_A,c);
            f_Value = f_Success(X_FinalPoint,Y_FinalPoint,X_Tx,Y_Tx,X_Rx,Y_Rx,Success_Set,SuccessCounter,Times_From_A,c);
            stepX = f_Value/dfx;
            stepY = f_Value/dfy;
            X_FinalPoint = X_FinalPoint - stepX;
            Y_FinalPoint = Y_FinalPoint - stepY;
            if(abs(X_FinalPoint)>MapLength || abs(Y_FinalPoint)>MapLength)
                SuccessSign = false;
                break;
            end
            
            ErrorX = (X_FinalPoint-Old_X_FinalPoint)/X_FinalPoint;
            ErrorY = (Y_FinalPoint-Old_Y_FinalPoint)/Y_FinalPoint;
            if(ErrorX^2 + ErrorY^2 < Error)
                X_FinalResult(time) = X_FinalPoint;
                Y_FinalResult(time) = Y_FinalPoint;
                SuccessSign = true;
                break;
            end
        end
        if(SuccessSign)
            break;
        end
    end
end
if(Time_Measure)
    t = toc;
    fprintf('Diff Method Time = %i\n',t);
end
end


%% Aux Functions
function [Out_dfx,Out_dfy] = find_dfx_dfy(X_A,Y_A,X_Tx,Y_Tx,X_Rx,Y_Rx,Success_Set,SuccessCounter,Times_From_A,c)
Out_dfx = 0;
Out_dfy = 0;
for NoRx = 1:SuccessCounter
    Out_dfx = Out_dfx + dfx(X_A,Y_A,X_Tx,Y_Tx,X_Rx,Y_Rx,NoRx,Success_Set,Times_From_A,c);
    Out_dfy = Out_dfy + dfy(X_A,Y_A,X_Tx,Y_Tx,X_Rx,Y_Rx,NoRx,Success_Set,Times_From_A,c);
end
Out_dfx = Out_dfx/SuccessCounter;
Out_dfy = Out_dfy/SuccessCounter;
end

function Out = f_Success(X_A,Y_A,X_Tx,Y_Tx,X_Rx,Y_Rx,Success_Set,SuccessCounter,Times_From_A,c)
Out = 0;
for NoRx = 1:SuccessCounter
    Out = Out + f(X_A,Y_A,X_Tx,Y_Tx,X_Rx,Y_Rx,NoRx,Success_Set,Times_From_A,c);
end
Out = Out/SuccessCounter;
end

function Out = f(X_A,Y_A,X_Tx,Y_Tx,X_Rx,Y_Rx,NoRx,Success_Set,Times_From_A,c)
Time_Temp = (sqrt((X_A-X_Rx(Success_Set(NoRx)))^2+(Y_A-Y_Rx(Success_Set(NoRx)))^2) + sqrt((X_A-X_Tx)^2+(Y_A-Y_Tx)^2))/c;
Out = (Time_Temp*10^7 - Times_From_A(Success_Set(NoRx))*10^7)^2;
end

function Out = dfx(X_A,Y_A,X_Tx,Y_Tx,X_Rx,Y_Rx,NoRx,Success_Set,Times_From_A,c)
delta = 0.00000001;
Out = (f(X_A + delta,Y_A,X_Tx,Y_Tx,X_Rx,Y_Rx,NoRx,Success_Set,Times_From_A,c) - f(X_A,Y_A,X_Tx,Y_Tx,X_Rx,Y_Rx,NoRx,Success_Set,Times_From_A,c))/delta;
end

function Out = dfy(X_A,Y_A,X_Tx,Y_Tx,X_Rx,Y_Rx,NoRx,Success_Set,Times_From_A,c)
delta = 0.00000001;
Out = (f(X_A,Y_A + delta,X_Tx,Y_Tx,X_Rx,Y_Rx,NoRx,Success_Set,Times_From_A,c) - f(X_A,Y_A,X_Tx,Y_Tx,X_Rx,Y_Rx,NoRx,Success_Set,Times_From_A,c))/delta;
end
