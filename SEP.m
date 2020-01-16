clc
clear all
xaxis=200;      %dimensionss of sensor network
yaxis=200;
m=0.2; %Percentage of nodes that are advance nodes
b=0.3;  %Percentage of nodes that are intermediate nodes

a=3;  %energy factor for advance nodes
u=a/2;  % energy factor for intermediate nodes

h=100;    %value for hard threshold
tEmpi=50;  
tEmpf=200;

Residual_e3(1)=0;
INFINITY = 999999999999999;

base_station.x=xaxis*0.5;  %distance of base station from the network
base_station.y=yaxis*0.5;
n = 200;         %no of nodes
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
maximum_lifetime=4000;           %no of rounds
do=sqrt(Energy_Frequency/Emplifier);       %distance between cluster head and base station

for i=1:1:n
    S(i).xd=rand(1,1)*xaxis;         %it will distribute the nodes in 1 dimension in x axis randomly.
    S(i).yd=rand(1,1)*yaxis;           %it will distribute the nodes in 1 dimension in y axis randoml.
    S(i).G=0;                        % initially there are no cluster heads
    S(i).E=initial_energy;
    S(i).type='N';                   %initially there are no cluster heads only nodes
end
S(n+1).xd=base_station.x;            %total no of nodes is n and with base station  it is n+1
S(n+1).yd=base_station.y;
countCHs=0;       
cluster=1;                            %first cluster is selected
flag_first_dead=0;
flag_teenth_dead=0;
flag_all_dead=0;
dead=0;                              % initially no node is dead
first_dead=0;
teenth_dead=0;
all_dead=0;                          
allive=n;                            %initially all nodes are alive
%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;
% counter for sleep nodes
s=0;

%Creation of the random Sensor Network

for i=1:1:n
    Se(i).xd=rand(1,1)*xaxis;
    XR(i)=S(i).xd;
    Se(i).yd=rand(1,1)*yaxis;
    YR(i)=S(i).yd;
    Se(i).G=0;
    %initially there are no cluster heads only nodes
    Se(i).type='N';
  
    temp_rnd0 = i;
    %Random Election of Normal Nodes
    if (temp_rnd0>=m*n+1)
        Se(i).E=initial_energy;
        Se(i).ENERGY=0;
       %%%plot(St(i).xd,St(i).yd,'o');
       hold on;
   end
    %Random Election of advance Nodes
   if  (temp_rnd0<m*n+1)
        Se(i).E=initial_energy*(1+a);
        Se(i).ENERGY=1;
       %%%%plot(St(i).xd,St(i).yd,'+');
         hold on;
   end
end
Se(n+1).xd=base_station.x;   %total no of nodes is n and with base station  it is n+1
Se(n+1).yd=base_station.y;
%plot(S2(n+1).xd,S2(n+1).yd,'x');
   
%counter for CHs
countCHs2=0;
%counter for CHs per round
rcountCHs2=0;
cluster2=1;

countCHs;
rcountCHs2=rcountCHs2+countCHs2;
flag_first_dead2=0;
alive=n;

for r=0:1:maximum_lifetime
    r
    warning('off','all');
  
  %Election Probability for Normal Nodes
    pnrm=( p/ (1+a*m) );
  %Election Probability for Advanced Nodes
 padv= ( p*(1+a)/(1+a*m) )
  
 %Operation for epoch 
 if(mod(r, round(1/pnrm) )==0)
    for i=1:1:n
        if(Se(i).ENERGY==0)
        Se(i).G=0;
        Se(i).cl=0;
    end
    end
  end
    
  
%Operation for epoch
if(mod(r, round(1/padv) )==0)
    for i=1:1:n
         if(Se(i).ENERGY==1)
        Se(i).G=0;
        Se(i).cl=0;
    end
    end
end


%Number of dead nodes
dead2=0;
%Number of dead Advanced Nodes
dead_a2=0;
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
    if (Se(i).E<=0)
        %plot(Se(i).xd,Se(i).yd,'red .');
        dead2=dead2+1;
        if(Se(i).ENERGY==1)
            dead_a2=dead_a2+1;
        end
        if(Se(i).ENERGY==0)
            dead_n2=dead_n2+1;
        % plot(Se(i).xd,Se(i).yd,'red .');
        end
            
     end
    if Se(i).E>0
        Se(i).type='N';
         if (Se(i).ENERGY==0)  
        %plot(Se(i).xd,Se(i).yd,'o');
         end
        if (Se(i).ENERGY==1)  
        %plot(Se(i).xd,Se(i).yd,'+');
        end
    
    end
end
%plot(Se(n+1).xd,Se(n+1).yd,'x');


STATISTICS.DEAD2(r+1)=dead2;
STATISTICS.ALIVE2(r+1)=alive-dead2;
DEAD2(r+1)=dead2;
DEAD_N2(r+1)=dead_n2;
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
   if(Se(i).E>0)
     temp_rand=rand;    
     if ( (Se(i).G)<=0)
        %Election of Cluster Heads
       
       if( ( Se(i).ENERGY==0 && ( temp_rand <= ( pnrm / ( 1 - pnrm * mod(r,round(1/pnrm)) )) ) )  )  % for normal nodes
            countCHs2=countCHs2+1;
            packets_TO_BS2=packets_TO_BS2+1;
            PACKETS_TO_BS2(r+1)=packets_TO_BS2;
            
            Se(i).type = 'C';
            Se(i).G = round(1/p)-1;
            C(cluster2).xd = Se(i).xd;
            C(cluster2).yd = Se(i).yd;
            %plot(Se(i).xd,Se(i).yd,'k*');
           
            distance=sqrt( (Se(i).xd-(Se(n+1).xd) )^2 + (Se(i).yd-(Se(n+1).yd) )^2 );
            C(cluster2).distance = distance;
            C(cluster2).id = i;
            X(cluster2)=Se(i).xd;
            Y(cluster2)=Se(i).yd;
            cluster2=cluster2+1;
           
            distanceBroad = sqrt(xaxis*xaxis+yaxis*yaxis);
            if (distanceBroad > do)
                Se(i).E = Se(i).E- ( (Energy_Transmitter)*(100) + Emplifier* 100*( distanceBroad*distanceBroad*distanceBroad*distanceBroad ));
            
            else
                Se(i).E = Se(i).E- ( (Energy_Transmitter)*(100) + Energy_Frequency * 100*( distanceBroad*distanceBroad));
            end
            %Calculation of Energy dissipated 
            
            distance;
           
            if (distance > do)
                 Se(i).E = Se(i).E- ( (Energy_Transmitter+Data_Aggregation_Energy)*(4000) + Emplifier * 4000*( distance*distance*distance*distance ));
            else
                 Se(i).E = Se(i).E- ( (Energy_Transmitter+Data_Aggregation_Energy)*(4000) + Energy_Frequency * 4000*( distance * distance ));
            end
        
            packets_TO_BS2 = packets_TO_BS2+1;
            PACKETS_TO_BS2(r+1) = packets_TO_BS2;
       end    
        
        if( ( Se(i).ENERGY==1 && ( temp_rand <= ( padv / ( 1 - padv * mod(r,round(1/padv)) )) ) )  )  % for advance nodes
            countCHs2=countCHs2+1;
            packets_TO_BS2=packets_TO_BS2+1;
            PACKETS_TO_BS2(r+1)=packets_TO_BS2;
            
            Se(i).type = 'C';
            Se(i).G = round(1/p)-1;
            C(cluster2).xd = Se(i).xd;
            C(cluster2).yd = Se(i).yd;
            %plot(Se(i).xd,Se(i).yd,'k*');
           
            distance=sqrt( (Se(i).xd-(Se(n+1).xd) )^2 + (Se(i).yd-(Se(n+1).yd) )^2 );
            C(cluster2).distance = distance;
            C(cluster2).id = i;
            X(cluster2)=Se(i).xd;
            Y(cluster2)=Se(i).yd;
            cluster2=cluster2+1;
           
            distanceBroad = sqrt(xaxis*xaxis+yaxis*yaxis);
            if (distanceBroad > do)
                Se(i).E = Se(i).E- ( (Energy_Transmitter)*(100) + Emplifier* 100*( distanceBroad*distanceBroad*distanceBroad*distanceBroad ));
            
            else
                Se(i).E = Se(i).E- ( (Energy_Transmitter)*(100) + Energy_Frequency * 100*( distanceBroad*distanceBroad));
            end
            %Calculation of Energy dissipated 
            distance;
          
            if (distance > do)
                 Se(i).E = Se(i).E- ( (Energy_Transmitter+Data_Aggregation_Energy)*(4000) + Emplifier * 4000*( distance*distance*distance*distance ));
            else
                 Se(i).E = Se(i).E- ( (Energy_Transmitter+Data_Aggregation_Energy)*(4000) + Energy_Frequency * 4000*( distance * distance ));
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
   if ( Se(i).type=='N' && Se(i).E>0 )
     if(cluster2-1>=1)
       min_dis=sqrt( (Se(i).xd-Se(n+1).xd)^2 + (Se(i).yd-Se(n+1).yd)^2 );
       min_dis_cluster2=1;
       for c=1:1:cluster2-1
           temp=min(min_dis,sqrt( (Se(i).xd-C(c).xd)^2 + (Se(i).yd-C(c).yd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster2=c;
           end
       end
       
       %Energy dissipated by associated Cluster Head
     
            min_dis;
            if (min_dis>do)
                Se(i).E=Se(i).E- ( Energy_Transmitter*(4000) + Emplifier*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                Se(i).E=Se(i).E- ( Energy_Transmitter*(4000) + Energy_Frequency*4000*( min_dis * min_dis)); 
            end
        %Energy dissipated
        if(min_dis>0)
            Se(C(min_dis_cluster2).id).E = Se(C(min_dis_cluster2).id).E- ( (Energy_Receiver + Data_Aggregation_Energy)*4000 ); 
         PACKETS_TO_CH2(r+1)=n-dead2-cluster2+1; 
        end
       
       Se(i).min_dis=min_dis;
       Se(i).min_dis_cluster2=min_dis_cluster2;
           
   end
 end
end

eee=0;
    for i=1:n
    eee=eee+Se(i).E;

    end
    
STATISTICS.PACKETS_TO_CH2(r+1)=packets_TO_CH2;
STATISTICS.PACKETS_TO_BS2(r+1)=packets_TO_BS2; 

    Residual_e3(r+1)=eee;         % residual energy of nodes

end