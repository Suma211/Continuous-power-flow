function [dif_PQ]=difference_PQ(P_sch,Q_sch,P_cal,Q_cal,PQ,nPQ) %difference between the scheduled and calculated values of PQ nodes.
dif_P_all=P_sch-P_cal;
dif_Q_all=Q_sch-Q_cal;

dif_Q=zeros(nPQ,1);
for i=1:nPQ
    dif_Q(i)=dif_Q_all(PQ(i));
end
    dif_P=dif_P_all(2:end);
    dif_PQ=[dif_P;dif_Q];
end
