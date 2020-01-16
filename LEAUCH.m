clc;
clear all;
close all;
%Sensor field X-Axis co-ordinates
xaxis=200;
%Sensor field Y-Axis co-ordinates
yaxis=200;
%Sink X-Axis co-ordinates
sink_x=250;
%Sink Y-Axis co-ordinates
sink_y=100;
%Number of sensor nodes
n=200;
%Probability to become cluster head
p1=0.04;
%Number of primary users
pu=5;
%Number of channels
k=5;
%Initial energy of the nodes
Eo=0.5;
%Energy consumption by the electronic model
Eelec=50e-9;
%Energy consumption by the free space model
Efs=10e-12;
%Energy consuption by the amplifier
Eamp=0.0013e-12;
%Energy consuption by data aggregation
Eda=5e-9;
%Distance
do=sqrt(Efs/Eamp);
%Maximum number of rounds
r=2000;
%Control message in 20 byte length
control_message=160;
%Data packet
data_packet=2000;
%Stores X-co-ordinates 
x1=zeros(n,1);
%Stores Y-co-ordinates
y1=zeros(n,1);
%Remaining energy of the network
Erem1=zeros(n,1);
%Stores distance between node and the sink
distance1=zeros(n,1);
%To store death of first node
first_dead1=0;
%To store when half of the network dies
half_dead1=0;
%To store when the last node dies
all_dead1=0;
%To break the loop when last node dies
flag1=0;
%Variable to count the number of dead nodes
dead_nodes1=0;
%Variable to count number of alive nodes
alive_nodes1=n;
%Variable to store packets sent to base station
packets_BS1=0;
%Variable to count the number of packet sent to cluster head
packets_CH1=0;
%Variable to count total number of packets
packet_count1=0;
%Variable to count the numebr of cluster heads formed in a round
CH_count1=0;
%Creation of Random sensor network field
figure(1);
for i=1:1:n
    %rand function generates a single random number between 0 & 1
    x1(i)=rand(1,1)*xaxis;
    y1(i)=rand(1,1)*yaxis;
    %plot is used to visualise the node deployment
    plot(x1(i),y1(i),'bo');
    hold on;
    %Setting the initial energy to all the nodes
    Erem1(i)=Eo;
    %Broadcast message by the sink to all the nodes
    Erem1(i)=Erem1(i)-(Eelec*160);
    %EPOCH
    S(i).G=0;
    %Distance between a node and the sink is found using Euclidean formula
    Distance1=sqrt((x1(i)-sink_x)^2+(y1(i)-sink_y)^2);
end
title('LEAUCH');
%Plot of the sink away from the field
plot(sink_x,sink_y,'hg','Markersize',20);
hold on;
for R1=1:1:r
    %Number of nodes that are alive after deployment
    alive_nodes1=n;
    for i=1:1:n
        %Condition to determine a dead node
        if(Erem1(i)<=0)
            %Increment dead node count by 1
            dead_node1=dead_node1+1;
            %Decrement alive node count by 1
            alive_node1=alive_node1-1;
        end
    end
    %Graph variable to plot remaining energy
    Residual_energy1(R1)=sum(Erem1);
    %Graph variable to plot total number of packets
    packets_count1(R1)=(packets_BS1+packets_CH1);
    %If flag variable is raised i.e 1 break the loop
    if(flag1==1)
        break;
    end
    %Variable to count number of cluster head
    cluster1=0;
    %Stores X-co-ordinates of cluster head
    Cx1=0;
    %Stores Y-co-ordinates of cluster head
    Cy1=0;
    %Stores distance to sink of a cluster head
    Cd1=0;
    %Stores ID of cluster head
    C_id1=0;
    %Stores remaining energy of cluster head
    Ce1=0;
    


