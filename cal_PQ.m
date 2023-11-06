
function [P_cal,Q_cal]=cal_PQ(V_mag,Y_mag,Theta,V_Delta,No_of_Buses) %calculate nodes'power flow

P_cal=zeros(No_of_Buses,1);  %initialize active power vector
Q_cal=zeros(No_of_Buses,1);  %initialize reactive power vector
     for i=1:No_of_Buses
         for j=1:No_of_Buses
             P_cal(i)=P_cal(i)+V_mag(i)*V_mag(j)*Y_mag(i,j)*cos(Theta(i,j)-V_Delta(i)+V_Delta(j));
             Q_cal(i)=Q_cal(i)-V_mag(i)*V_mag(j)*Y_mag(i,j)*sin(Theta(i,j)-V_Delta(i)+V_Delta(j));
         end
     end
end


