function [Y_mat,Theta,Y_mag,B,G]=y_bus(Bus_data,Line_data,No_of_Buses,No_of_Lines)
Y_mat=zeros(No_of_Buses);
B_mat=zeros(No_of_Buses);
Z_x=Line_data(:,2); %one side bus number of branches/transformers 
Z_y=Line_data(:,3); %another side bus number of branches/transformers
R=Line_data(:,8);  %value of R
X=Line_data(:,9);  %value of X
B=Line_data(:,10).*1j; %value of B
t_trans=Line_data(:,16); %transformer final turns ratio
Bus_shunt=Bus_data(:,15)+Bus_data(:,16)*1i;
Z=R+1j*X; %impedance
Y=1./Z; %admittance
for n = 1:No_of_Lines    
    if t_trans(n)~= 0    %off-diagnoal elements of transformer admittance matrix   
        Y_mat(Z_x(n),Z_y(n))=-1/t_trans(n)*Y(n);
        Y_mat(Z_y(n),Z_x(n))=Y_mat(Z_x(n), Z_y(n));
        Y_mat(Z_x(n),Z_x(n))=Y_mat(Z_x(n),Z_x(n))+Y(n)*(1/t_trans(n))^2;
        Y_mat(Z_y(n),Z_y(n))=Y_mat(Z_y(n),Z_y(n))+Y(n);         
    else
        Y_mat(Z_x(n),Z_y(n))=-Y(n); %off-diagnoal elements of branch admittance matrix
        Y_mat(Z_y(n),Z_x(n))=Y_mat(Z_x(n),Z_y(n));  
        Y_mat(Z_x(n),Z_x(n))=Y_mat(Z_x(n),Z_x(n))+Y(n)+B(n)/2;
        Y_mat(Z_y(n),Z_y(n))=Y_mat(Z_y(n),Z_y(n))+Y(n)+B(n)/2;   
    end
end
for k=1:No_of_Buses %diagnoal elements plus bus admittance
    Y_mat(k,k)=Y_mat(k,k)+Bus_shunt(k);
end
[Theta, Y_mag] = cart2pol(real(Y_mat), imag(Y_mat)); % get angle and magnitude of admittance matrix
 B = imag(Y_mat); % susceptance matrix
 G = real(Y_mat); % conductance matrix
end




