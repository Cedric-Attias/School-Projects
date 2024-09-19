% Project: KIN 613 - 3D Kinematics Assignment
% Cedric Attias
%--------------------------------------------------------------------------
% Perpare MATLAB to run script by closing open windsows and clearing both
% the workspace and command window

close all
clear
clc

%--------------------------------------------------------------------------
%% PREPARING DATA
%--------------------------------------------------------------------------
% Crop the Data Manually Outside of Matlab
% Remove columns (Physical marker data (Cols 2-49)
% Remove empty rows at the beggining (Varies by trial)
% Remove empty rows at the end (Varies by trial)

% Import CSV file into MATLAB
%FilePathGaitSample = "ExampleGaitData_Cropped.csv";
FilePathGait = "Gait.csv";
FilePath_Squat = "Squat.csv";
FilePath_Stand = "Stand.csv";
FilePath_UniKneel = "UnilateralKneeling.csv";
FilePath_SupKneel = "SingleArmSupportedKneeling.csv";

dataArray = readmatrix(FilePath_SupKneel);

% Define Constants
freq = 64; % sampling frequncy (Hz)
deltaT = 1/freq; % time increment based on sampling frequency
[rowsSize, colsSize] = size(dataArray); %Number of rows and cols

%Cubic spline interpolation
interpolatedData = zeros(rowsSize, colsSize);
interpolatedData(:,1)= dataArray(:,1);
for j = 2:colsSize
    interpolatedData(:,j) = fillmissing(dataArray(:,j),'spline');
end

%Check interpolated data against Raw Data of LFCON Y-Trajectory (Quick
%Check to amke sure data is continuous --> Needed for differentiation)
%figure
%plot(interpolatedData(:,1), interpolatedData(:,24),'b', dataArray(:,1), dataArray(:,24), 'r', 'linewidth', 5);
%title('Check Interpolation: LFCON Y-Trajectory ');

dataArray = interpolatedData;

% Filter data set using a 2nd-order lowpass Butterworth filter with cutoff
% frequency of 6 Hz (Litertaure in Biomechanical work for simple movements)

cf = 6; % cut-off frequency

[b,a] = butter(2,cf/(freq*0.5),"low");

for i = 2:colsSize
    dataArray(:,i) = filtfilt(b,a,dataArray(:,i)); 
end

%Check filtered data against Raw Data of LFCON Y-Trajectory (Because data has holes)
%figure
%plot(interpolatedData(:,1), interpolatedData(:,24),'b', dataArray(:,1), dataArray(:,24), 'r', 'linewidth', 5);
%title('Check Interpolation: LFCON Y-Trajectory ');
%--------------------------------------------------------------------------
%% DATA ANALYSIS
%--------------------------------------------------------------------------
% Define key variables

frames = 1:size(dataArray,1); %Frame array used to loop through rows

timeSeries = (frames-1)*deltaT; %Time series for position
timeSeriesDer =  (frames(2:end-1)-1)*deltaT; %Time series for derived values (vel, acc)

% For Global Coordinate System (GCS): ----
% Define unit vectors

I = [1 0 0]; % Global X-axis unit vector
J = [0 1 0]; % Global Y-axis unit vector
K = [0 0 1]; % Global Z-axis unit vector

% Global to local rotation matrix (ZYX roation)

%Local to Global ZYX rotation matrix --> ISB Standard

%L_R_G = [(cosd(alpha)*cosd(beta)) ((cosd(alpha)*sind(beta)*sind(gamma))-(sind(alpha)*cosd(gamma))) ((cosd(alpha)*sind(beta)*cosd(gamma))+(sind(alpha)*sind(gamma)));
%(sind(alpha)*cosd(beta)) ((sind(alpha)*sind(beta)*sind(gamma))+(cosd(alpha)*cosd(gamma))) ((sind(alpha)*sind(beta)*cosd(gamma))-(cosd(alpha)*sind(gamma)));
%(-sind(beta)) (cosd(beta)*sind(gamma)) (cosd(beta)*cosd(gamma))];

%beta = -asind(dot(K,i));
%alpha = (asind(dot(J,i)))/(cosd(beta));
%gamma = (asind(dot(K,j)))/(cosd(beta));

% Separate data into arrays for each landamrk 

P_R_ASIS    = dataArray(:,2:4);
P_R_ILIAC_C = dataArray(:,5:7);
P_R_PSIS    = dataArray(:,8:10);
P_L_PSIS    = dataArray(:,11:13);
P_L_ILIAC_C = dataArray(:,14:16);
P_L_ASIS    = dataArray(:,17:19);
T_R_GT      = dataArray(:,20:22);
T_R_L_FCON  = dataArray(:,23:25);
T_R_M_FCON  = dataArray(:,26:28);
S_R_M_TCON  = dataArray(:,29:31);
S_R_L_TCON  = dataArray(:,32:34);
S_R_TIB_TUB = dataArray(:,35:37);
S_R_M_MAL   = dataArray(:,38:40);
S_R_L_MAL   = dataArray(:,41:43);
F_R_L_MAL   = dataArray(:,44:46);
F_R_M_MAL   = dataArray(:,47:49);
F_R_H       = dataArray(:,50:52);
F_R_T       = dataArray(:,53:55);
F_R_B5_MT   = dataArray(:,56:58);
F_R_H1_MT   = dataArray(:,59:61);
%--------------------------------------------------------------------------
% Calculate ankle joint center, AJC

AJC = zeros(length(frames),3); %Variable pre-allocation!!

for j=1:length(frames)
    AJC(j,:) = (F_R_L_MAL(j,:) + F_R_M_MAL(j,:))/2;
end
%--------------------------------------------------------------------------
% Calculate knee joint center, KJC

KJC = zeros(length(frames),3);

for j=1:length(frames)
    KJC(j,:) = (T_R_L_FCON(j,:) + T_R_M_FCON(j,:))/2;
end
%--------------------------------------------------------------------------
% Calculate hip joint center, HJC

% initialize variables
Width_P = zeros(length(frames),1);
O_P = zeros(length(frames),3);
Mid_PSIS = zeros(length(frames),3);

z_P = zeros(length(frames),3);
k_P = zeros(length(frames),3);

temp_v_P = zeros(length(frames),3);

y_P = zeros(length(frames),3);
j_P = zeros(length(frames),3);

x_P = zeros(length(frames),3);
i_P = zeros(length(frames),3);

beta_P = zeros(length(frames),1);
alpha_P = zeros(length(frames),1);
gamma_P = zeros(length(frames),1);

HJC_global = zeros(length(frames),3);


for j=1:length(frames)   
    %Pelvis width from R to L ASIS
    Width_P(j) = norm(P_R_ASIS(j,:) - P_L_ASIS(j,:));
    
    % Calculate origin of Pelvis in the LCS
    O_P(j,:) = (P_L_ASIS(j,:) + P_R_ASIS(j,:))/2;
     
    %Find the midpoint of the L and R PSIS
    Mid_PSIS(j,:) = (P_L_PSIS(j,:) + P_R_PSIS(j,:))/2;
    
    % Define local axes for pevlis LCS and local unit vectors
    % z-axis:
    z_P(j,:) = P_R_ASIS(j,:) - O_P(j,:);
    k_P(j,:) = z_P(j,:)/norm(z_P(j,:)); %Unit vector
    
    %temporary axis:
    temp_v_P(j,:) = O_P(j,:) - Mid_PSIS(j,:);
    
    % y-axis:
    y_P(j,:) = cross(z_P(j,:),temp_v_P(j,:));
    j_P(j,:) = y_P(j,:)/norm(y_P(j,:));%Unit vector
    
    %x-axis
    x_P(j,:) = cross(y_P(j,:),z_P(j,:));
    i_P(j,:) = x_P(j,:)/norm(x_P(j,:));
    
    % Calculate gamma, beta, and alpha for rotation matrix
    % beta:
    beta_P(j,:) = -asind(dot(K,i_P(j,:)));

    % alpha:
    alpha_P(j,:) = asind(dot(J,i_P(j,:))/cosd(beta_P(j,:)));

    % gamma
    gamma_P(j,:) = asind(dot(K,j_P(j,:))/cosd(beta_P(j,:)));
     
end

%Take average of pelvis width coordinates
Width_P_Avg = mean(Width_P);

%The HJC in the LCS (Based on literature)
HJC_local = [(-0.19*Width_P_Avg), (-0.30*Width_P_Avg), (0.36*Width_P_Avg)];

% Calculate HJC in GCS

for j=1:length(frames)
   
   %Pelvis to Global Rotation Matrix
   P_R_G = [(cosd(alpha_P(j))*cosd(beta_P(j))) ((cosd(alpha_P(j))*sind(beta_P(j))*sind(gamma_P(j)))-(sind(alpha_P(j))*cosd(gamma_P(j)))) ((cosd(alpha_P(j))*sind(beta_P(j))*cosd(gamma_P(j)))+(sind(alpha_P(j))*sind(gamma_P(j))));
   (sind(alpha_P(j))*cosd(beta_P(j))) ((sind(alpha_P(j))*sind(beta_P(j))*sind(gamma_P(j)))+(cosd(alpha_P(j))*cosd(gamma_P(j)))) ((sind(alpha_P(j))*sind(beta_P(j))*cosd(gamma_P(j)))-(cosd(alpha_P(j))*sind(gamma_P(j))));
   (-sind(beta_P(j))) (cosd(beta_P(j))*sind(gamma_P(j))) (cosd(beta_P(j))*cosd(gamma_P(j)))];

   HJC_global(j,:) = O_P(j,:) + transpose((P_R_G)*transpose(HJC_local));
end
%--------------------------------------------------------------------------
% Calculate velocity and acceleration for joint centers
% Using 3-point central difference method to calculate both
% Derivation would cause additional noise

v_AJC = zeros(length(frames)-2,3);
v_KJC = zeros(length(frames)-2,3);
v_HJC = zeros(length(frames)-2,3);
a_AJC = zeros(length(frames)-2,3);
a_KJC = zeros(length(frames)-2,3);
a_HJC = zeros(length(frames)-2,3);

for j=2:length(frames)-1
    % Velocities
    v_AJC(j-1,:) = (AJC(j+1,:) - AJC(j-1,:))/(2*deltaT);
    v_KJC(j-1,:) = (KJC(j+1,:) - KJC(j-1,:))/(2*deltaT);
    v_HJC(j-1,:) = (HJC_global(j+1,:) - HJC_global(j-1,:))/(2*deltaT);
    % Accelerations
    a_AJC(j-1,:) = (AJC(j+1,:) - 2*AJC(j,:) + AJC(j-1,:))/deltaT^2;
    a_KJC(j-1,:) = (KJC(j+1,:) - 2*KJC(j,:) + KJC(j-1,:))/deltaT^2;
    a_HJC(j-1,:) = (HJC_global(j+1,:) - 2*HJC_global(j,:) + HJC_global(j-1,:))/deltaT^2;
end

%--------------------------------------------------------------------------
% Calculate LCS for foot
O_F = zeros(length(frames),3); %Foot origin (Heel)

z_F = zeros(length(frames),3);
k_F = zeros(length(frames),3);

temp_v_F = zeros(length(frames),3);

y_F = zeros(length(frames),3);
j_F = zeros(length(frames),3);

x_F = zeros(length(frames),3);
i_F = zeros(length(frames),3);

Mid_Mal_F = zeros(length(frames),3);

for j=1:length(frames)
    O_F(j,:) = F_R_H(j,:); %Foot origin (Heel) (ISB standard)
    
    % Define local axes for the foot LCS and local unit vectors
    % x-axis:Line going through the 2nd metatarsal from the heel origin 
    x_F(j,:) = F_R_T(j,:) - O_F(j,:);
    i_F(j,:) = x_F(j,:)/norm(x_F(j,:)); %Unit vector
    
    %Find the midpoint of the Malleoli
    Mid_Mal_F(j,:) = (F_R_M_MAL(j,:) + F_R_L_MAL(j,:))/2;
    
    %temporary axis: Heel to mid mal
    temp_v_F(j,:) = Mid_Mal_F(j,:) - O_F(j,:); %Pointing backwards towards the heel
    
    % z-axis:
    z_F(j,:) = cross(x_F(j,:),temp_v_F(j,:));
    k_F(j,:) = z_F(j,:)/norm(z_F(j,:));%Unit vector
    
    %y-axis
    y_F(j,:) = cross(z_F(j,:),x_F(j,:));
    j_F(j,:) = y_F(j,:)/norm(y_F(j,:)); 
end
%--------------------------------------------------------------------------
% Calculate LCS for shank
z_S = zeros(length(frames),3);
k_S = zeros(length(frames),3);

temp_v_S = zeros(length(frames),3);

y_S = zeros(length(frames),3);
j_S = zeros(length(frames),3);

x_S = zeros(length(frames),3);
i_S = zeros(length(frames),3);

MID_Mal_S = ones(length(frames),3);

O_S = zeros(length(frames),3); 

for j = 1:length(frames)
    
    O_S(j,:) = ((T_R_M_FCON(j,:))+(T_R_L_FCON(j,:)))/2; %Origin @ midpoint of femoral condyles 
    
    MID_Mal_S(j,:) = (S_R_M_MAL(j,:) + S_R_L_MAL(j,:))/2;
    
    % Define local axes for the shank LCS and local unit vectors
    % y-axis: Line through the origin of the shank from the malleoli midpoint
    y_S(j,:) = O_S(j,:) -  MID_Mal_S(j,:); 
    j_S(j,:) = y_S(j,:)/norm(y_S(j,:)); %Unit vector pointing up the longitudical axis of the shank
    
    %temporary axis: between medial to lateral fem cond
    temp_v_S(j,:) = T_R_L_FCON(j,:) - T_R_M_FCON(j,:);
    
    % x-axis:
    x_S(j,:) = cross(y_S(j,:),temp_v_S(j,:));
    i_S(j,:) = x_S(j,:)/norm(x_S(j,:));%Unit vector
    
    % z-axis:
    z_S(j,:) = cross(x_S(j,:),y_S(j,:));
    k_S(j,:) = z_S(j,:)/norm(z_S(j,:));%Unit vector
end
%--------------------------------------------------------------------------
% Calculate LCS for Thigh

z_T = zeros(length(frames),3);
k_T = zeros(length(frames),3);

temp_v_T = zeros(length(frames),3);

y_T = zeros(length(frames),3);
j_T = zeros(length(frames),3);

x_T = zeros(length(frames),3);
i_T = zeros(length(frames),3);

MID_FCON = ones(length(frames),3);

O_T = zeros(length(frames),3); 

for j = 1:length(frames)
    
    O_T(j,:) = HJC_global(j,:); %Origin of the thigh is at the hip joint centre (ISB standard)
    
    MID_FCON(j,:) = ((T_R_M_FCON(j,:))+(T_R_L_FCON(j,:)))/2; %Midpoint of Fem Conds
   
    % Define local axes for the thigh LCS and local unit vectors
    % y-axis: Line through the origin of the thigh from the mid fem conds
    y_T(j,:) = O_T(j,:) -  MID_FCON(j,:); 
    j_T(j,:) = y_T(j,:)/norm(y_T(j,:)); %Unit vector
    
    %temporary axis: between medial to lateral fem cond
    temp_v_T(j,:) = T_R_L_FCON(j,:) - T_R_M_FCON(j,:);
    
    % x-axis:
    x_T(j,:) = cross(y_T(j,:),temp_v_T(j,:));
    i_T(j,:) = x_T(j,:)/norm(x_T(j,:));%Unit vector
    
    % z-axis:
    z_T(j,:) = cross(x_T(j,:),y_T(j,:));
    k_T(j,:) = z_T(j,:)/norm(z_T(j,:));%Unit vector
end
%--------------------------------------------------------------------------
% Calculate Euler angles for foot, shank, and thigh for local to global
% rotation of each segment

alpha_F = ones(length(frames),1);
beta_F = ones(length(frames),1);
gamma_F = ones(length(frames),1);

alpha_S = ones(length(frames),1);
beta_S = ones(length(frames),1);
gamma_S = ones(length(frames),1);

alpha_T = ones(length(frames),1);
beta_T = ones(length(frames),1);
gamma_T = ones(length(frames),1);

for j=1:length(frames)    
    % For foot:
     beta_F(j,:) = -asind(dot(K,i_F(j,:)));
     alpha_F(j,:) = asind(dot(J,i_F(j,:))/cosd(beta_F(j,:)));
     gamma_F(j,:)   = asind(dot(K,j_F(j,:))/cosd(beta_F(j,:)));

    % For shank:
    beta_S(j,:) = -asind(dot(K,i_S(j,:)));
    alpha_S(j,:) = asind(dot(J,i_S(j,:))/cosd(beta_S(j,:)));
    gamma_S(j,:)   = asind(dot(K,j_S(j,:))/cosd(beta_S(j,:)));

    % For thigh:
    beta_T(j,:) = -asind(dot(K,i_T(j,:)));
    alpha_T(j,:) = asind(dot(J,i_T(j,:))/cosd(beta_T(j,:)));
    gamma_T(j,:)   = asind(dot(K,j_T(j,:))/cosd(beta_T(j,:)));
end

%--------------------------------------------------------------------------
% Intersegmental angles for ankle (between foot and shank = SF) and knee (between shank and thigh = TS)
alpha_TS = ones(length(frames),1);
gamma_TS = ones(length(frames),1);
beta_TS = ones(length(frames),1);
e2_TS = ones(length(frames),3);

alpha_SF = ones(length(frames),1);
gamma_SF = ones(length(frames),1);
beta_SF = ones(length(frames),1);
e2_SF = ones(length(frames),3);

for j=1:length(frames)   
    
    % For knee joint:
    % Floating axis method versus using an additional rotation matrix
    e2_TS(j,:) = cross(j_S(j,:),k_T(j,:))/norm(cross(j_S(j,:),k_T(j,:)));
    alpha_TS(j,:) = asind(dot(-1*e2_TS(j,:),j_T(j,:))); %flexion (+)/extension(-) - Z axis
    beta_TS(j,:) = asind(dot(-1*e2_TS(j,:),k_S(j,:))); % iternal(+)/external(-) rotation - Y axis
    gamma_TS(j,:) = acosd(dot(k_T(j,:),j_S(j,:))) - 90; %adduction(+)/abduction(-) - X axis

    % For ankle joint:
    e2_SF(j,:) = cross(j_S(j,:),k_F(j,:))/norm(cross(j_S(j,:),k_F(j,:)));
    alpha_SF(j,:) = asind(dot(-1*e2_SF(j,:),j_F(j,:))); %dorsiflexion(+)/plantarflexion(-) - Z
    beta_SF(j,:) = asind(dot(-1*e2_SF(j,:),k_S(j,:))); %adduction(+)/abduction(-) - Y
    gamma_SF(j,:) = acosd(dot(k_F(j,:),j_S(j,:))) - 90; %inversion(+)/eversion(-) - X
end

%--------------------------------------------------------------------------
%% RESULTS 
%--------------------------------------------------------------------------
% For static trials - find averages and standard deviations of angles

% Segment angles:
% Foot:
alpha_F_mean = mean(alpha_F);
alpha_F_sd = std(alpha_F);
beta_F_mean = mean(beta_F);
beta_F_sd = std(beta_F);
gamma_F_mean = mean(gamma_F);
gamma_F_sd = std(gamma_F);

% Shank:
alpha_S_mean = mean(alpha_S);
alpha_S_sd = std(alpha_S);
beta_S_mean = mean(beta_S);
beta_S_sd = std(beta_S);
gamma_S_mean = mean(gamma_S);
gamma_S_sd = std(gamma_S);

% Thigh:
alpha_T_mean = mean(alpha_T);
alpha_T_sd = std(alpha_T);
beta_T_mean = mean(beta_T);
beta_T_sd = std(beta_T);
gamma_T_mean = mean(gamma_T);
gamma_T_sd = std(gamma_T);

% Joint angles:

% Ankle:
alpha_SF_mean = mean(alpha_SF);
alpha_SF_sd = std(alpha_SF);
beta_SF_mean = mean(beta_SF);
beta_SF_sd = std(beta_SF);
gamma_SF_mean = mean(gamma_SF);
gamma_SF_sd = std(gamma_SF);

% Knee:
alpha_TS_mean = mean(alpha_TS);
alpha_TS_sd = std(alpha_TS);
beta_TS_mean = mean(beta_TS);
beta_TS_sd = std(beta_TS);
gamma_TS_mean = mean(gamma_TS);
gamma_TS_sd = std(gamma_TS);

% Put calculations into array/table for easy interpretation and
% visualization

SegmentAngles = ["Segment Angle", "Angle Mean", "Angle SD";
                 "Foot A",  alpha_F_mean, alpha_F_sd;
                 "Foot B",  beta_F_mean,  beta_F_sd;
                 "Foot G",  gamma_F_mean,   gamma_F_sd;
                 "Shank A", alpha_S_mean, alpha_S_sd;
                 "Shank B", beta_S_mean,  beta_S_sd;
                 "Shank G", gamma_S_mean,   gamma_S_sd;
                 "Thigh A", alpha_T_mean, alpha_T_sd;
                 "Thigh B", beta_T_mean,  beta_T_sd;
                 "Thigh G", gamma_T_mean,   gamma_T_sd;];

JointAngles = ["Joint Angle", "Angle Mean", "Angle SD";
               "Ankle A", alpha_SF_mean, alpha_SF_sd;
               "Ankle B", beta_SF_mean,  beta_SF_sd;
               "Ankle G", gamma_SF_mean,   gamma_SF_sd;
               "Knee A",  alpha_TS_mean, alpha_TS_sd;
               "Knee B",  beta_TS_mean,  beta_TS_sd;
               "Knee G",  gamma_TS_mean,   gamma_TS_sd;];
%--------------------------------------------------------------------------
% For walking trial - plot data

% Plot position, velocity and acceleration for ankle joint centre (all directions)

figure
subplot(3,3,1)
plot(timeSeries,AJC(:,1));
title('AJC X-Displacement in GCS')
xlabel('Time (s)')
ylabel('X Position (mm)')
grid on
subplot(3,3,2)
plot(timeSeries,AJC(:,2));
title('AJC Y-Displacement in GCS')
xlabel('Time (s)')
ylabel('Y Position (mm)')
grid on
subplot(3,3,3)
plot(timeSeries,AJC(:,3));
title('AJC Z-Displacement in GCS')
xlabel('Time (s)')
ylabel('Z Position (mm)')
grid on

subplot(3,3,4)
plot(timeSeriesDer,v_AJC(:,1));
title('AJC X-Velocity in GCS')
xlabel('Time (s)')
ylabel('X Velocity (mm/s)')
grid on
subplot(3,3,5)
plot(timeSeriesDer,v_AJC(:,2));
title('AJC Y-Velocity in GCS')
xlabel('Time (s)')
ylabel('Y Velocity (mm/s)')
grid on
subplot(3,3,6)
plot(timeSeriesDer,v_AJC(:,3));
title('AJC Z-Velocity in GCS')
xlabel('Time (s)')
ylabel('Z Velocity (mm/s)')
grid on

subplot(3,3,7)
plot(timeSeriesDer,a_AJC(:,1));
title('AJC X-Acceleration in GCS')
xlabel('Time (s)')
ylabel('X Acceleration (mm/s^2)')
grid on
subplot(3,3,8)
plot(timeSeriesDer,a_AJC(:,2));
title('AJC Y-Acceleration in GCS')
xlabel('Time (s)')
ylabel('Y Acceleration (mm/s^2)')
grid on
subplot(3,3,9)
plot(timeSeriesDer,a_AJC(:,3));
title('AJC Z-Acceleration in GCS')
xlabel('Time (s)')
ylabel('Z Acceleration GCS (mm/s^2)')
grid on
%--------------------------------------------------------------------------
% Plot position, velocity and acceleration for knee joint centre (all directions)

figure
subplot(3,3,1)
plot(timeSeries,KJC(:,1));
title('KJC X-Displacement in GCS')
xlabel('Time (s)')
ylabel('X Position (mm)')
grid on
subplot(3,3,2)
plot(timeSeries,KJC(:,2));
title('KJC Y-Displacement in GCS')
xlabel('Time (s)')
ylabel('Y Position (mm)')
grid on
subplot(3,3,3)
plot(timeSeries,KJC(:,3));
title('KJC Z-Displacement in GCS')
xlabel('Time (s)')
ylabel('Z Position (mm)')
grid on

subplot(3,3,4)
plot(timeSeriesDer,v_KJC(:,1));
title('KJC X-Velocity in GCS')
xlabel('Time (s)')
ylabel('X Velocity (mm/s)')
grid on
subplot(3,3,5)
plot(timeSeriesDer,v_KJC(:,2));
title('KJC Y-Velocity in GCS')
xlabel('Time (s)')
ylabel('Y Velocity (mm/s)')
grid on
subplot(3,3,6)
plot(timeSeriesDer,v_KJC(:,3));
title('KJC Z-Velocity in GCS')
xlabel('Time (s)')
ylabel('Z Velocity (mm/s)')
grid on

subplot(3,3,7)
plot(timeSeriesDer,a_KJC(:,1));
title('KJC X-Acceleration in GCS')
xlabel('Time (s)')
ylabel('X Acceleration (mm/s^2)')
grid on
subplot(3,3,8)
plot(timeSeriesDer,a_KJC(:,2));
title('KJC Y-Acceleration in GCS')
xlabel('Time (s)')
ylabel('Y Acceleration (mm/s^2)')
grid on
subplot(3,3,9)
plot(timeSeriesDer,a_KJC(:,3));
title('KJC Z-Acceleration in GCS')
xlabel('Time (s)')
ylabel('Z Acceleration (mm/s^2)')
grid on
%--------------------------------------------------------------------------
% Plot position, velocity and acceleration for hip joint centre (all directions)

figure
subplot(3,3,1)
plot(timeSeries,HJC_global(:,1))
title('HJC X-Displacement in GCS')
xlabel('Time (s)')
ylabel('X Position (mm)')
grid on
subplot(3,3,2)
plot(timeSeries,HJC_global(:,2))
title('HJC Y-Displacement in GCS')
xlabel('Time (s)')
ylabel('Y Position (mm)')
grid on
subplot(3,3,3)
plot(timeSeries,HJC_global(:,3))
title('HJC Z-Displacement in GCS')
xlabel('Time (s)')
ylabel('Z Position (mm)')
grid on

subplot(3,3,4)
plot(timeSeriesDer,v_HJC(:,1));
title('HJC X-Velocity in GCS')
xlabel('Time (s)')
ylabel('X Velocity (mm/s)')
grid on
subplot(3,3,5)
plot(timeSeriesDer,v_HJC(:,2));
title('HJC Y-Velocity in GCS')
xlabel('Time (s)')
ylabel('Y Velocity (mm/s)')
grid on
subplot(3,3,6)
plot(timeSeriesDer,v_HJC(:,3));
title('HJC Z-Velocity in GCS')
xlabel('Time (s)')
ylabel('Z Velocity (mm/s)')
grid on

subplot(3,3,7)
plot(timeSeriesDer,a_HJC(:,1));
title('HJC X-Acceleration in GCS')
xlabel('Time (s)')
ylabel('X Acceleration (mm/s^2)')
grid on
subplot(3,3,8)
plot(timeSeriesDer,a_HJC(:,2));
title('HJC Y-Acceleration in GCS')
xlabel('Time (s)')
ylabel('Y Acceleration (mm/s^2)')
grid on
subplot(3,3,9)
plot(timeSeriesDer,a_HJC(:,3));
title('HJC Z -Acceleration in GCS')
xlabel('Time (s)')
ylabel('Z Acceleration (mm/s^2)')
grid on
%--------------------------------------------------------------------------
% Plot walking segment angle data

figure
subplot(3,3,1)
plot(timeSeries,alpha_T);
title('Thigh Alpha Angle w.r.t. GCS')
xlabel('Time (s)')
ylabel('Thigh Alpha (deg)')
grid on
subplot(3,3,2)
plot(timeSeries,beta_T);
title('Thigh Beta Angle w.r.t. GCS')
xlabel('Time (s)')
ylabel('Thigh Beta (deg)')
grid on
subplot(3,3,3)
plot(timeSeries,gamma_T);
title('Thigh Gamma Angle w.r.t. GCS')
xlabel('Time (s)')
ylabel('Thigh Gamma (deg)')
grid on

subplot(3,3,4)
plot(timeSeries,alpha_S);
title('Shank Alpha Angle w.r.t. GCS')
xlabel('Time (s)')
ylabel('Shank Alpha (deg)')
grid on
subplot(3,3,5)
plot(timeSeries,beta_S);
title('Shank Beta Angle w.r.t. GCS')
xlabel('Time (s)')
ylabel('Shank Beta (deg)')
grid on
subplot(3,3,6)
plot(timeSeries,gamma_S);
title('Shank Gamma Angle w.r.t. GCS')
xlabel('Time (s)')
ylabel('Shank Gamma (deg)')
grid on

subplot(3,3,7)
plot(timeSeries,alpha_F);
title('Foot Alpha Angle w.r.t. GCS')
xlabel('Time (s)')
ylabel('Foot Alpha (deg)')
grid on
subplot(3,3,8)
plot(timeSeries,beta_F);
title('Foot Beta Angle w.r.t. GCS')
xlabel('Time (s)')
ylabel('Foot Beta (deg)')
grid on
subplot(3,3,9)
plot(timeSeries,gamma_F);
title('Foot Gamma Angle w.r.t. GCS')
xlabel('Time (s)')
ylabel('Foot Gamma (deg)')
grid on
%--------------------------------------------------------------------------
% Segment angles - Ankle and knee

figure
subplot(1,3,1)
plot(timeSeries,alpha_SF);
title('Ankle Dorsi/Plantarflexion (Alpha)')
xlabel('Time (s)')
ylabel('Ankle PF(-)/DF(+) (deg)')
grid minor
subplot(1,3,2)
plot(timeSeries,gamma_SF);
title('Ankle Inversion/Eversion (Gamma)')
xlabel('Time (s)')
ylabel('Ankle EV(-)/INV(+) (deg)')
grid minor
subplot(1,3,3)
plot(timeSeries,beta_SF);
title('Ankle Abduction/Adduction (Beta)')
xlabel('Time (s)')
ylabel('Ankle ABD(-)/ADD(+) (deg)')
grid minor

figure
subplot(1,3,1)
plot(timeSeries,alpha_TS);
title('Knee Flexion/Extension (Alpha)')
xlabel('Time (s)')
ylabel('Knee EXT(-)/FLEX(+) (deg)')
grid minor
subplot(1,3,2)
plot(timeSeries,gamma_TS);
title('Knee Abduction/Adduction (Gamma)')
xlabel('Time (s)')
ylabel('Knee ABD(-)/ADD(+) (deg)')
grid minor
subplot(1,3,3)
plot(timeSeries,beta_TS);
title('Knee Internal/External Rotation (Beta)')
xlabel('Time (s)')
ylabel('Knee EXT(-)/INT(+) (deg)')
grid minor