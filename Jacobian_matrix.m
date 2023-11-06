%% Jacobian matrix
function [J]=Jacobian_matrix(V_mag,P_cal,Q_cal,Y_mag,Theta, V_Delta,No_of_Buses,PQ,nPQ,B,G)

J1=zeros(No_of_Buses-1);
J2=zeros(No_of_Buses-1,nPQ);

J3=zeros(nPQ,No_of_Buses-1);
J4=zeros(nPQ,nPQ);


%J1=dP/dDelta
for i=2:No_of_Buses %row position of a bus
    for j=2:No_of_Buses %column position of bus
        if(j==i) 
                       
            J1(i-1,j-1)=-Q_cal(i)-V_mag(i)^2*B(i,j); %diagonal elements  
        else
            J1(i-1,j-1)= -V_mag(i)*V_mag(j)*Y_mag(i,j)*sin(Theta(i,j)-V_Delta(i)+V_Delta(j));;
        end
    end
end

 %J2=dP/dV
 for i=2:No_of_Buses %position of row
     for j=1:nPQ  %position of column
         jj=PQ(j);         
         if(jj==i) %diagonal elements     
             J2(i-1,j)=P_cal(i)/V_mag(i)+V_mag(i)*G(i,i);
            
         else %off-diagonal elements
            J2(i-1,j)= V_mag(i)*Y_mag(i,jj)*cos(Theta(i,jj)-V_Delta(i)+V_Delta(jj));
                 
         end
     end
 end
             
 %J3=dQ/dDelta
 for i=1:nPQ %position of row
     ii=PQ(i);
     for j=2:No_of_Buses %position of column
         if(j==ii) 
             J3(i,j-1)=P_cal(ii)-G(ii,ii)*(V_mag(ii)^2); %diagonal elements
            
         else  %off-diagonal elements
             J3(i,j-1)=-V_mag(ii)*V_mag(j)*Y_mag(ii,j)*cos(Theta(ii,j)-V_Delta(ii)+V_Delta(j));
         end
     end
 end
 
 %J4=dQ/dV
 for i=1:nPQ %position of row
     ii=PQ(i);
     for j=1:nPQ %position of colume
         jj=PQ(j);
         if(jj==ii) %diagonal elements 
             J4(i,j)=Q_cal(ii)/V_mag(ii)-V_mag(ii)*B(ii,ii);
             
         else %off-diagonal elements
             J4(i,j)=-V_mag(ii)*Y_mag(ii,jj)*sin(Theta(ii,jj)-V_Delta(ii)+V_Delta(jj));
                
         end
     end
 end
   J=[J1 J2;J3 J4];  
 end
 

