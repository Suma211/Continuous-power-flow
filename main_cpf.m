clc; 
clear all
%% read grid data from file
file_name='ieee14cdf.txt';
[S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data]=read_data(file_name);

 V_mag_final=Bus_data(:,5); %magnitude of final voltage in p.u.
 V_ang_final=Bus_data(:,6)*pi/180; % change angle of final voltage from degrees to radians.
 P_load=Bus_data(:,7)/S_Base; %active power of Load in p.u.
 Q_load=Bus_data(:,8)/S_Base; %reactive power of Load in p.u.
 P_gen=Bus_data(:,9)/S_Base; %active power of generation in p.u.
 Q_gen=Bus_data(:,10)/S_Base; %reactive power of genertion in p.u.
 Qmax = Bus_data(:,13)/S_Base;     % Maximum Reactive Power Limit
 Qmin = Bus_data(:,14)/S_Base;      % Minimum Reactive Power Limit

 %% type of buses
 slack=find(Bus_data(:,4)==3); %slack bus
 PQ=find(Bus_data(:,4)==0|Bus_data (:,4) == 1); %PQ bus index
 PV=find(Bus_data(:,4)==2); %PV bus index
 nslack=length(slack);
 nPQ=length(PQ);
 nPV=length(PV);
 %% form Y_matrix
 [Y_mat,Theta,Y_mag,B,G]=y_bus(Bus_data,Line_data,No_of_Buses,No_of_Lines); %form Y_matrix
 
%% initialize parameters
P_sch=P_gen-P_load; %net power scheduled at a bus
Q_sch=Q_gen-Q_load;
V_mag = V_mag_final;  
V_Delta = V_ang_final;   
Lambda=0;
CPF_prediction=[V_Delta(2:end);V_mag(PQ); Lambda];
CPF_correction=[V_Delta(2:end);V_mag(PQ); Lambda]; 

%% phase 1 upper part of PV curve
%predictor
while 1
    sigma=0.1;
    K=[P_sch(2:end);Q_sch(PQ)];
    [P_cal,Q_cal]=cal_PQ(V_mag,Y_mag,Theta,V_Delta,No_of_Buses); %calculate the power at buses     
    [J]=Jacobian_matrix(V_mag,P_cal,Q_cal,Y_mag,Theta, V_Delta,No_of_Buses,PQ,nPQ,B,G); %call the Jacobian matrix
   ek=zeros(1,23);
   ek(end)=1;
   J_extended=[J,-K;ek];
   predict_vector=[V_Delta(2:end);V_mag(PQ); Lambda]+sigma*inv(J_extended)*ek';  
   V_Delta(2:end)=predict_vector(1:nPQ+nPV);
   V_mag(PQ)=predict_vector(nPQ+nPV+1:end-1);
   Lambda=predict_vector(end);
   CPF_prediction=[CPF_prediction predict_vector];  

 %corrector--solve power flow equations by NR method

 tol_max=0.01; % iteration tolerance of the solution
 iter=0; %times of iteration
 tol=1; %initial value of tol
  P_sch=(1+Lambda)*(P_gen-P_load); %load change-since Lambda will change
  Q_sch=(1+Lambda)*(Q_gen-Q_load);
 while (tol>tol_max && iter<10)    
    [P_cal,Q_cal]=cal_PQ(V_mag,Y_mag,Theta,V_Delta,No_of_Buses); %calculate the power at buses     
    [dif_PQ]=difference_PQ(P_sch,Q_sch,P_cal,Q_cal,PQ,nPQ); % mismatches vector
    [J]=Jacobian_matrix(V_mag,P_cal,Q_cal,Y_mag,Theta, V_Delta,No_of_Buses,PQ,nPQ,B,G); %call the Jacobian matrix
    dif_Voltage=inv(J)*dif_PQ; %get correction vector
    dif_D=dif_Voltage(1:No_of_Buses-1); % angle correction vector
    dif_V=dif_Voltage(No_of_Buses:end); % magnitude correction vector
    
    V_Delta(2:end)=V_Delta(2:end)+dif_D; %correct the results, angle and voltage
    V_mag(PQ)=V_mag(PQ)+dif_V;
    tol=max(abs(dif_PQ));
    iter=iter+1;
 end
correct_vector=[V_Delta(2:end);V_mag(PQ); Lambda];
CPF_correction =[CPF_correction correct_vector];
 if iter>=10
       CPF_correction(:,end)=[]; % if not iteration, the value of last column got from non-converged calculation 
       CPF_prediction(:,end)=[]; 
       break;
 end 
 end
CPF_correction(:,end)=[]; 
CPF_prediction(:,end)=[];
V_Delta(2:end)=CPF_correction(1:13,end);
V_mag(PQ)=CPF_correction(14:22,end);
Lambda=CPF_correction(end,end);
P_sch = (1+Lambda)*(P_gen-P_load);  
Q_sch = (1+Lambda)*(Q_gen-Q_load);  
       
figure;
plot(CPF_prediction(22,:),'-r');
hold on
plot(CPF_correction(22,:),'-b');

legend('V after predictor','V after corrector')

figure;
plot(CPF_correction(end,:),CPF_correction(14:22,:),'.-')

%% phase 2  near the tip of the nose curve--switch from changing lamda to V
 counter=0;
while 1
    sigma=0.025;
    K=[P_sch(2:end);Q_sch(PQ)];    
    [P_cal,Q_cal]=cal_PQ(V_mag,Y_mag,Theta,V_Delta,No_of_Buses); %calculate the power at buses     
    [J]=Jacobian_matrix(V_mag,P_cal,Q_cal,Y_mag,Theta, V_Delta,No_of_Buses,PQ,nPQ,B,G); %call the Jacobian matrix
   ek=zeros(1,23);
   ek(1,23)=1;
   ekv=zeros(1,23);
   ekv(1,22)=-1;
   J_extended=[J,-K;ekv]; 
   
   predict_vector=[V_Delta(2:end);V_mag(PQ); Lambda]+sigma*inv(J_extended)*ek';  
   V_Delta(2:end)=predict_vector(1:nPQ+nPV);
   V_mag(PQ)=predict_vector(nPQ+nPV+1:end-1);
   Lambda=predict_vector(end);
   CPF_prediction=[CPF_prediction predict_vector];  

%% corrector--solve power flow equations by NR method
tol_max=0.01; % iteration tolerance of the solution
 iter=0; %times of iteration
 tol=1; %initial value of tol
 
 while (tol>tol_max && iter<10)
     P_sch=(1+Lambda)*(P_gen-P_load); %load change
     Q_sch=(1+Lambda)*(Q_gen-Q_load);
    [P_cal,Q_cal]=cal_PQ(V_mag,Y_mag,Theta,V_Delta,No_of_Buses); %calculate the power at buses     
    [dif_PQ]=difference_PQ(P_sch,Q_sch,P_cal,Q_cal,PQ,nPQ); % mismatches vector
    [J]=Jacobian_matrix(V_mag,P_cal,Q_cal,Y_mag,Theta, V_Delta,No_of_Buses,PQ,nPQ,B,G); %call the Jacobian matrix
    
    J_extended=[J,K;ek];
    M=[dif_PQ;0];
    X=inv(J_extended)*M; 
    
    dif_D=X(1:No_of_Buses-1); %angle correction vector
    dif_V=X(No_of_Buses:end-1); %magnitude correction vector
    
    V_Delta(2:end)=V_Delta(2:end)+dif_D; %correct the results, angle and voltage
    V_mag(PQ)=V_mag(PQ)+dif_V;
    Lambda=Lambda+X(end);
    
    tol=max(abs(dif_PQ));
    iter=iter+1;
 end
correct_vector=[V_Delta(2:end);V_mag(PQ); Lambda];
CPF_correction =[CPF_correction correct_vector];

 if iter>=10
       CPF_correction(:,end)=[]; % if not iteration, the value of last column got from non-converged calculation 
       CPF_prediction(:,end)=[]; 
       break;
 end 

 if CPF_correction(end,end)<CPF_correction(end,end-1)
    counter=counter+1;
 end
    if counter ==5
        break;
    end
end
V_Delta(2:end)=CPF_correction(1:13,end);
V_mag(PQ)=CPF_correction(14:22,end);
Lambda=CPF_correction(end,end);
P_sch = (1+Lambda)*(P_gen-P_load);  
Q_sch = (1+Lambda)*(Q_gen-Q_load);  

%% phase 3 the lower part of PV curve
% predictor
while 1
    sigma=0.1;
    K=[P_sch(2:end);Q_sch(PQ)];
    [P_cal,Q_cal]=cal_PQ(V_mag,Y_mag,Theta,V_Delta,No_of_Buses); %calculate the power at buses     
    [J]=Jacobian_matrix(V_mag,P_cal,Q_cal,Y_mag,Theta, V_Delta,No_of_Buses,PQ,nPQ,B,G); %call the Jacobian matrix
   ek=zeros(1,23);
ek=zeros(1,23);
ek(1,23)=1;
ekv=zeros(1,23);
ekv(1,23)=-1;
  
   J_extended=[J,-K;ekv];
   
   predict_vector=[V_Delta(2:end);V_mag(PQ); Lambda]+sigma*inv(J_extended)*ek';  
   V_Delta(2:end)=predict_vector(1:nPQ+nPV);
   V_mag(PQ)=predict_vector(nPQ+nPV+1:end-1);
   Lambda=predict_vector(end);
   CPF_prediction=[CPF_prediction predict_vector];  

 %corrector--solve power flow equations by NR method
P_sch=(1+Lambda)*(P_gen-P_load); %load change
Q_sch=(1+Lambda)*(Q_gen-Q_load);

tol_max=0.01; % iteration tolerance of the solution
 iter=0; %times of iteration
 tol=1; %initial value of tol
 
 while (tol>tol_max && iter<10)
    [P_cal,Q_cal]=cal_PQ(V_mag,Y_mag,Theta,V_Delta,No_of_Buses); %calculate the power at buses     
    [dif_PQ]=difference_PQ(P_sch,Q_sch,P_cal,Q_cal,PQ,nPQ); % mismatches vector
    [J]=Jacobian_matrix(V_mag,P_cal,Q_cal,Y_mag,Theta, V_Delta,No_of_Buses,PQ,nPQ,B,G); %call the Jacobian matrix
    dif_Voltage=inv(J)*dif_PQ; %get correction vector
    dif_D=dif_Voltage(1:No_of_Buses-1); % angle correction vector
    dif_V=dif_Voltage(No_of_Buses:end); % magnitude correction vector
    
    V_Delta(2:end)=V_Delta(2:end)+dif_D; %correct the results, angle and voltage
    V_mag(PQ)=V_mag(PQ)+dif_V;
    tol=max(abs(dif_PQ));
    iter=iter+1;
 end
correct_vector=[V_Delta(2:end);V_mag(PQ); Lambda];
CPF_correction =[CPF_correction correct_vector];
 
  if CPF_correction(end,end)<=0
       break;
   end   
end 
figure;
plot(CPF_correction(end,:),CPF_correction(14,:),'.-');
xlabel('Lambda');
ylabel('Voltage magnitude of bus 4 (pu)');

figure;
plot(CPF_correction(end,:),CPF_correction(15,:),'.-');
xlabel('Lambda');
ylabel('Voltage magnitude of bus 5 (pu)');

figure;
plot(CPF_correction(end,:),CPF_correction(16,:),'.-');
xlabel('Lambda');
ylabel('Voltage magnitude of bus 7 (pu)');

figure;
plot(CPF_correction(end,:),CPF_correction(17,:),'.-');
xlabel('Lambda');
ylabel('Voltage magnitude of bus 9 (pu)');

figure;
plot(CPF_correction(end,:),CPF_correction(18,:),'.-');
xlabel('Lambda');
ylabel('Voltage magnitude of bus 10 (pu)');

figure;
plot(CPF_correction(end,:),CPF_correction(19,:),'.-');
xlabel('Lambda');
ylabel('Voltage magnitude of bus 11 (pu)');

figure;
plot(CPF_correction(end,:),CPF_correction(20,:),'.-');
xlabel('Lambda');
ylabel('Voltage magnitude of bus 12 (pu)');

figure;
plot(CPF_correction(end,:),CPF_correction(21,:),'.-');
xlabel('Lambda');
ylabel('Voltage magnitude of bus 13 (pu)');

figure;
plot(CPF_correction(end,:),CPF_correction(22,:),'.-');
xlabel('Lambda');
ylabel('Voltage magnitude of bus 14 (pu)');

