clc
close all
clear all;
xaxis=200;      %dimensionss of sensor network
yaxis=200;
base_station.x=0.5*xaxis;  %distance of base station from the network
base_station.y=1.5*yaxis;
n = 100;         %no of nodes
p= 0.05;          %probibilty of a node to become cluster head
initial_energy=0.5;          %energy supplied to each node
Energy_Transmitter=50*0.000000001;     %energy dissipated per bit at transmitter per node
Energy_Receiver=50*0.000000001;        %energy dissipated per bit at reciever per node
Energy_Frequency=10*0.000000000001;    %amplification coefficient of free-space signal when d is less than d0
Emplifier=0.0013*0.000000000001;      % multi-path fading signal amplification coefficient when d is greater than d0
Residual_e1(1)=0;            
Residual_e2(1)=0;
Residual_e3(1)=0;
Residual_e4(1)=0;
Residual_e5(1)=0;
Data_Aggregation_Energy=5*0.000000001;           % data aggregation energy per node
maximum_lifetime=9999;           %no of rounds
do=sqrt(Energy_Frequency/Emplifier);       %distance between cluster head and base station
for i=1:1:n
    S(i).xd=rand(1,1)*xaxis;         %it will distribute the nodes in 1 dimension in x axis randomly.
    S(i).yd=rand(1,1)*yaxis;           %it will distribute the nodes in 1 dimension in y axis randoml.
    S(i).G=0;                        % initially there are no cluster heads
    S(i).E=initial_energy;
    S(i).type='N';                   %initially there are no cluster heads only nodes
end
%Values for Hetereogeneity

m=0.2; %Percentage of nodes that are advance nodes
b=0.3;  %Percentage of nodes that are intermediate nodes

a=3;  %energy factor for advance nodes
u=a/2;  % energy factor for intermediate nodes

h=100;    %value for hard threshold
tEmpi=50;  
tEmpf=200;

Residual_e3(1)=0;
INFINITY = 999999999999999;

%Creation of the random Sensor Network

for i=1:1:n
    St(i).xd=rand(1,1)*xaxis;
    XR(i)=S(i).xd;
    St(i).yd=rand(1,1)*yaxis;
    YR(i)=S(i).yd;
    St(i).G=0;
    %initially there are no cluster heads only nodes
    S(i).type='N';
  
    temp_rnd0 = i;
    %Random Election of Normal Nodes
    if (temp_rnd0>=(b+m)*n+1)
        St(i).E=initial_energy;
        St(i).ENERGY=0;
       %%%plot(St(i).xd,St(i).yd,'o');
       hold on;
   end
    %Random Election of intermediate Nodes
   if (temp_rnd0<(b+m)*n+1)&& (temp_rnd0>=m*n+1)
        St(i).E=initial_energy*(1+u);
        St(i).ENERGY=0.5;
       %%%%plot(St(i).xd,St(i).yd,'+');
         hold on;
   end
    %Random Election of Advanced Nodes
   if (temp_rnd0<m*n+1)
        St(i).E=initial_energy*(1+a);
        St(i).ENERGY=1;
        hold on;
   end
end
St(n+1).xd=base_station.x;
St(n+1).yd=base_station.y;
%plot(St(n+1).xd,St(n+1).yd,'x');
   
%counter for CHs
countCHs2=0;
%counter for CHs per round
rcountCHs2=0;
cluster2=1;


rcountCHs2=rcountCHs2+countCHs2;
flag_first_dead2=0;
alive=n;

for r=0:1:maximum_lifetime
    r
    warning('off','all');
   
    cv = tEmpi + (tEmpf-tEmpi).*rand(1,1);
  %Election Probability for Normal Nodes
  pnrm=( p/ (1+a*m+u*b) );
  %Election Probability for intermediate Nodes
   pint=( p*(1+u)/ (1+a*m+u*b) );
  %Election Probability for Advanced Nodes
  padv= ( p*(1+a)/(1+a*m+u*b) );
 
  %Operation for epoch
  if(mod(r, round(1/pnrm) )==0)
    for i=1:1:n
        if(St(i).ENERGY==0)
        St(i).G=0;
        St(i).cl=0;
    end
    end
  end
    
    if(mod(r, round(1/pint) )==0)
    for i=1:1:n
        if(St(i).ENERGY==0.5)
        St(i).G=0;
        St(i).cl=0;
    end
    end
    end
  
if(mod(r, round(1/padv) )==0)
    for i=1:1:n
         if(St(i).ENERGY==1)
        St(i).G=0;
        St(i).cl=0;
    end
    end
end


%Number of dead nodes
dead2=0;
%Number of dead Advanced Nodes
dead_a2=0;
%number of intermediate nodes
dead_i2=0;
%Number of dead Normal Nodes
dead_n2=0;



%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS2=0;
packets_TO_CH2=0;
%counter for bit transmitted to Bases Station and to Cluster Heads
%per round
PACKETS_TO_CH2(r+1)=0;
PACKETS_TO_BS2(r+1)=0;


for i=1:1:n
    %checking if there is a dead node
    if (St(i).E<=0)
        %plot(St(i).xd,St(i).yd,'red .');
        dead2=dead2+1;
        if(St(i).ENERGY==1)
            dead_a2=dead_a2+1;
        end
        if(St(i).ENERGY==.5)
            dead_i2=dead_i2+1;
        end
        if(St(i).ENERGY==0)
            dead_n2=dead_n2+1;
        % plot(S(i).xd,S(i).yd,'red .');
        end
            
     end
    if St(i).E>0
        St(i).type='N';
         if (St(i).ENERGY==0)  
        %plot(S(i).xd,S(i).yd,'o');
         end
        if (St(i).ENERGY==0.5)  
        %plot(St(i).xd,St(i).yd,'^');
        end
        if (St(i).ENERGY==1)  
        %plot(S(i).xd,S(i).yd,'+');
        end
    
    end
end
%plot(S(n+1).xd,S(n+1).yd,'x');


STATISTICS.DEAD2(r+1)=dead2;
STATISTICS.ALIVE(r+1)=alive-dead2;
DEAD2(r+1)=dead2;
DEAD_N2(r+1)=dead_n2;
DEAD_I2(r+1)=dead_i2;
DEAD_A2(r+1)=dead_a2;



%When the first node dies
if (dead2==1)
    if(flag_first_dead2==0)
        first_dead2=r
        flag_first_dead2=1;
    end
end


countCHs2=0;
cluster2=1;
for i=1:1:n
   if(St(i).E>0)
     temp_rand=rand;    
     if ( (St(i).G)<=0)
        %Election of Cluster Heads
       
       if( ( St(i).ENERGY==0 && ( temp_rand <= ( pnrm / ( 1 - pnrm * mod(r,round(1/pnrm)) )) ) )  )  % for normal nodes
            countCHs2=countCHs2+1;
            packets_TO_BS2=packets_TO_BS2+1;
            PACKETS_TO_BS2(r+1)=packets_TO_BS2;
            
            St(i).type = 'C';
            St(i).G = round(1/p)-1;
            C(cluster2).xd = St(i).xd;
            C(cluster2).yd = St(i).yd;
            %plot(St(i).xd,St(i).yd,'k*');
           
            distance=sqrt( (St(i).xd-(St(n+1).xd) )^2 + (St(i).yd-(St(n+1).yd) )^2 );
            C(cluster2).distance = distance;
            C(cluster2).id = i;
            X(cluster2)=St(i).xd;
            Y(cluster2)=St(i).yd;
            cluster2=cluster2+1;
           
            distanceBroad = sqrt(xaxis*xaxis+yaxis*yaxis);
            if (distanceBroad > do)
                St(i).E = St(i).E- ( (Energy_Transmitter)*(100) + Emplifier* 100*( distanceBroad*distanceBroad*distanceBroad*distanceBroad ));
            
            else
                St(i).E = St(i).E- ( (Energy_Transmitter)*(100) + Energy_Frequency * 100*( distanceBroad*distanceBroad));
            end
            %Calculation of Energy dissipated 
            
            distance;
            if(cv>=h)
            if (distance > do)
                 St(i).E = St(i).E- ( (Energy_Transmitter+Data_Aggregation_Energy)*(4000) + Emplifier * 4000*( distance*distance*distance*distance ));
            else
                 St(i).E = St(i).E- ( (Energy_Transmitter+Data_Aggregation_Energy)*(4000) + Energy_Frequency * 4000*( distance * distance ));
            end
            end
            packets_TO_BS2 = packets_TO_BS2+1;
            PACKETS_TO_BS2(r+1) = packets_TO_BS2;
       end    
        if( ( St(i).ENERGY==0.5 && ( temp_rand <= ( pint / ( 1 - pint * mod(r,round(1/pint)) )) ) )  )  % for intermediate nodes
            countCHs2=countCHs2+1;
            packets_TO_BS2=packets_TO_BS2+1;
            PACKETS_TO_BS2(r+1)=packets_TO_BS2;
            
            St(i).type = 'C';
            St(i).G = round(1/p)-1;
            C(cluster2).xd = St(i).xd;
            C(cluster2).yd = St(i).yd;
            %plot(St(i).xd,St(i).yd,'k*');
           
            distance=sqrt( (St(i).xd-(St(n+1).xd) )^2 + (St(i).yd-(St(n+1).yd) )^2 );
            C(cluster2).distance = distance;
            C(cluster2).id = i;
            X(cluster2)=St(i).xd;
            Y(cluster2)=St(i).yd;
            cluster2=cluster2+1;
           
            distanceBroad = sqrt(xaxis*xaxis+yaxis*yaxis);
            if (distanceBroad > do)
                St(i).E = St(i).E- ( (Energy_Transmitter)*(100) + Emplifier* 100*( distanceBroad*distanceBroad*distanceBroad*distanceBroad ));
            
            else
                St(i).E = St(i).E- ( (Energy_Transmitter)*(100) + Energy_Frequency * 100*( distanceBroad*distanceBroad));
            end
            %Calculation of Energy dissipated 
            distance;
            if(cv>=h)      % if current sensed value is greater than hard threshold then only transmission will be done 
            if (distance > do)
                 St(i).E = St(i).E- ( (Energy_Transmitter+Data_Aggregation_Energy)*(4000) + Emplifier * 4000*( distance*distance*distance*distance ));
            else
                 St(i).E = St(i).E- ( (Energy_Transmitter+Data_Aggregation_Energy)*(4000) + Energy_Frequency * 4000*( distance * distance ));
            end
            end
            packets_TO_BS2 = packets_TO_BS2+1;
            PACKETS_TO_BS2(r+1) = packets_TO_BS2;
        end    
        if( ( St(i).ENERGY==1 && ( temp_rand <= ( padv / ( 1 - padv * mod(r,round(1/padv)) )) ) )  )  % for advance nodes
            countCHs2=countCHs2+1;
            packets_TO_BS2=packets_TO_BS2+1;
            PACKETS_TO_BS2(r+1)=packets_TO_BS2;
            
            St(i).type = 'C';
            St(i).G = round(1/p)-1;
            C(cluster2).xd = St(i).xd;
            C(cluster2).yd = St(i).yd;
            %plot(St(i).xd,St(i).yd,'k*');
           
            distance=sqrt( (St(i).xd-(St(n+1).xd) )^2 + (St(i).yd-(St(n+1).yd) )^2 );
            C(cluster2).distance = distance;
            C(cluster2).id = i;
            X(cluster2)=St(i).xd;
            Y(cluster2)=St(i).yd;
            cluster2=cluster2+1;
           
            distanceBroad = sqrt(xaxis*xaxis+yaxis*yaxis);
            if (distanceBroad > do)
                St(i).E = St(i).E- ( (Energy_Transmitter)*(100) + Emplifier* 100*( distanceBroad*distanceBroad*distanceBroad*distanceBroad ));
            
            else
                St(i).E = St(i).E- ( (Energy_Transmitter)*(100) + Energy_Frequency * 100*( distanceBroad*distanceBroad));
            end
            %Calculation of Energy dissipated 
            distance;
            if(cv>=h)
            if (distance > do)
                 St(i).E = St(i).E- ( (Energy_Transmitter+Data_Aggregation_Energy)*(4000) + Emplifier * 4000*( distance*distance*distance*distance ));
            else
                 St(i).E = St(i).E- ( (Energy_Transmitter+Data_Aggregation_Energy)*(4000) + Energy_Frequency * 4000*( distance * distance ));
            end
            end
            packets_TO_BS2 = packets_TO_BS2+1;
            PACKETS_TO_BS2(r+1) = packets_TO_BS2;
        end    
   
     end
   end
end

    STATISTICS.COUNTCHS2(r+1)=countCHs2;

STATISTICS.CLUSTERHEADS2(r+1) = cluster2-1;
CLUSTERHS2(r+1)= cluster2-1;

for i=1:1:n
   if ( St(i).type=='N' && St(i).E>0 )
     if(cluster2-1>=1)
       min_dis=sqrt( (St(i).xd-St(n+1).xd)^2 + (St(i).yd-St(n+1).yd)^2 );
       min_dis_cluster2=1;
       for c=1:1:cluster2-1
           temp=min(min_dis,sqrt( (St(i).xd-C(c).xd)^2 + (St(i).yd-C(c).yd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster2=c;
           end
       end
       
       %Energy dissipated by associated Cluster Head
       if(cv>=h)
            min_dis;
            if (min_dis>do)
                St(i).E=St(i).E- ( Energy_Transmitter*(4000) + Emplifier*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                St(i).E=St(i).E- ( Energy_Transmitter*(4000) + Energy_Frequency*4000*( min_dis * min_dis)); 
            end
        %Energy dissipated
        if(min_dis>0)
            St(C(min_dis_cluster2).id).E = St(C(min_dis_cluster2).id).E- ( (Energy_Receiver + Data_Aggregation_Energy)*4000 ); 
         PACKETS_TO_CH2(r+1)=n-dead2-cluster2+1; 
        end
       end
       St(i).min_dis=min_dis;
       St(i).min_dis_cluster2=min_dis_cluster2;
           
   end
 end
end

eee=0;
    for i=1:n
    eee=eee+St(i).E;

    end
    
STATISTICS.PACKETS_TO_CH2(r+1)=packets_TO_CH2;
STATISTICS.PACKETS_TO_BS2(r+1)=packets_TO_BS2; 

    Residual_e6(r+1)=eee;             % residual energy of nodes

end
