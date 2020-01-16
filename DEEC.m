%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEEC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all
clear all
xm=200;
ym=200;
sink.x=xm*0.5;  %location of sink on x-axis
sink.y=ym*0.5;  %location of sink on y-axis
n=200  %nodes
P=0.1;  %probability of cluster heads
Eo=0.5; %initial energy
ETX=50*0.000000001;  %tx energy
ERX=50*0.000000001;  %rx energy
Efs=10*0.000000000001;  %free space loss
Emp=0.0013*0.000000000001;   %multi path loss
%Data Aggregation Energy
EDA=5*0.000000001;  %compression energy
a=1;   %fraction of energy enhancment of advance nodes
rmax=4000  %maximum number of rounds
do=sqrt(Efs/Emp);  %distance do is measured
Et=0;  %variable just use below 
A=0;
for i=1:1:n
    S(i).xd=rand(1,1)*xm;  %generates a random no. use to randomly distibutes nodes on x axis
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;  %generates a random no. use to randomly distibutes nodes on y axis
    YR(i)=S(i).yd;
    S(i).G=0; %node is elegible to become cluster head
    talha=rand*a;
    S(i).E=Eo*(1+talha);
    E(i)= S(i).E;
    A=A+talha;
    Et=Et+E(i);  %estimating total energy of the network
    %initially there are no cluster heads only nodes
    S(i).type='N';
end
d1=0.765*xm/2;  %distance between cluster head and base station
K=sqrt(0.5*n*do/pi)*xm/d1^2; %optimal no. of cluster heads
d2=xm/sqrt(2*pi*K);  %distance between cluster members and cluster head
Er=4000*(2*n*ETX+n*EDA+K*Emp*d1^4+n*Efs*d2^2);  %energy desipated in a round
S(n+1).xd=sink.x; %sink is a n+1 node, x-axis postion of a node
S(n+1).yd=sink.y; %sink is a n+1 node, y-axis postion of a node
countCHs=0;  %variable, counts the cluster head
cluster=1;  %cluster is initialized as 1
flag_first_dead=0; %flag tells the first node dead
flag_teenth_dead=0;  %flag tells the 10th node dead
flag_all_dead=0;  %flag tells all nodes dead
dead=0;  %dead nodes count initialized to 0
first_dead=0;
teenth_dead=0;
all_dead=0;
allive=n;
%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;
for r=0:1:rmax     
    r
  if(mod(r, round(1/P) )==0)
    for i=1:1:n
        S(i).G=0;
        S(i).cl=0;
    end
  end
Ea=Et*(1-r/rmax)/n;
dead=0;
for i=1:1:n
   
    if (S(i).E<=0)
        dead=dead+1; 
        if (dead==1)
           if(flag_first_dead==0)
              first_dead=r;
              flag_first_dead=1;
           end
        end
        if(dead==0.1*n)
           if(flag_teenth_dead==0)
              teenth_dead=r;
              flag_teenth_dead=1;
           end
        end
        if(dead==n)
           if(flag_all_dead==0)
              all_dead=r;
              flag_all_dead=1;
           end
        end
    end
    if S(i).E>0
        S(i).type='N';
    end
end
STATISTICS.DEAD(r+1)=dead;
STATISTICS.ALLIVE(r+1)=allive-dead;
countCHs=0;
cluster=1;


for i=1:1:n
 if Ea>0
    
     
  p(i)=P*n*(1+a)*E(i)/(n+A)*(Ea);   
 %p(i)=P*n*S(i).E*E(i)/(Et*Ea);
 
 
 if(S(i).E>0)
   temp_rand=rand;     
   if ( (S(i).G)<=0)  
        if(temp_rand<= (p(i)/(1-p(i)*mod(r,round(1/p(i))))))
            countCHs=countCHs+1;
            packets_TO_BS=packets_TO_BS+1;
            PACKETS_TO_BS(r+1)=packets_TO_BS;
             S(i).type='C';
            S(i).G=round(1/p(i))-1;
            C(cluster).xd=S(i).xd;
            C(cluster).yd=S(i).yd;
           distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
            C(cluster).distance=distance;
            C(cluster).id=i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
            cluster=cluster+1;
           distance;
            if (distance>do)
                S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); 
            end
            if (distance<=do)
                S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance )); 
            end
        end     
    
   end
   
 end 
 end
end
STATISTICS.COUNTCHS(r+1)=countCHs;
%(5)簇内成员选择簇头模块(即簇的形成模块)
%簇内成员对簇头的选择（即簇的形成）算法
for i=1:1:n
   if ( S(i).type=='N' && S(i).E>0 )
     if(cluster-1>=1)
       min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
       min_dis_cluster=0;
       for c=1:1:cluster-1
           temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end
       %簇内节点（发送4000bit数据）能量消耗
       if(min_dis_cluster~=0)    
            min_dis;
            if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
      
            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
            packets_TO_CH=packets_TO_CH+1;
       else 
            min_dis;
            if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
            packets_TO_BS=packets_TO_BS+1;
            
       end
        S(i).min_dis=min_dis;
       S(i).min_dis_cluster=min_dis_cluster;
   else
            min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
            if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
            packets_TO_BS=packets_TO_BS+1;
   end
  end
end
STATISTICS.PACKETS_TO_CH(r+1)=packets_TO_CH;
STATISTICS.PACKETS_TO_BS(r+1)=packets_TO_BS;
RES(r+1)=0;
for i=1:1:n
    RES(r+1)=RES(r+1)+S(i).E;
end
end
first_dead
teenth_dead
all_dead
% STATISTICS.DEAD(r+1)
% STATISTICS.ALLIVE(r+1)
% STATISTICS.PACKETS_TO_CH(r+1)
% STATISTICS.PACKETS_TO_BS(r+1)
% STATISTICS.COUNTCHS(r+1)
% r=0:5000;
% subplot(2,2,1);
% plot(STATISTICS.DEAD,r);
% subplot(2,2,2);
% plot(r,STATISTICS.ALLIVE);
% subplot(2,2,3);
% plot(r,STATISTICS.PACKETS_TO_BS);
% subplot(2,2,4);
% plot(r,STATISTICS.COUNTCHS);

x=1:1:r;
y=1:1:r;
w=1:1:r;
v=1:1:r;

for i=1:r;
    x(i)=i;
    y(i) = n- STATISTICS.DEAD(i);   % remaining number of live nodes
    w(i)=RES(i);
    v(i)=STATISTICS.PACKETS_TO_BS(i);
end

figure(1)
plot(x,y,'r');
hold on;

figure(2)
plot(x,w,'r');
hold on;

figure(3)
plot(x,v,'r');
hold on;
