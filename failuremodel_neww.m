%  close all;
%  clear all;

 function failuremodel_new(topology,k_per)
 
% % n1=500;  %number of nodes in network X
% % n2=500;  %number of nodes in network Y
% p1=(4/(n1-1)); %the probability of having a link between nodes of network X
% p2=(4/(n2-1)); %the probability of having a link between nodes of network Y
% p12=0.05; %prob. of having a link from node X to node Y
% p21=0.05; %prob. of having a link from node Y to node X

% % variables for Barabasi Arbert Model/Scale Free Network
% m10 = 3; %m10: number of initially placed nodes of network X
% m1  = 2; %m11: number of nodes a new added node is connected to, 1 <= m1 < m10
% m20 = 3; %m20: number of initially placed nodes of network Y
% m2  = 2; %m21: number of nodes a new added node is connected to, 1 <= m2 < m20

% % variables for Watts Strogatz Model/Small Workd Network
% deg1 = 4; %d1: mean degree of network X, deg%2=0 && 0<deg<n-1
% rp1 = 0.5; %r1: rewiring probability of network X, 0<=rp<=1
% deg2 = 4; %d2: mean degree of network Y, deg%2=0 && 0<deg<n-1 
% rp2 = 0.5; %r2: rewiring probability of network Y, 0<=rp<=1


% %variables for Forest Fire model
% pFF=0.1; %forward burning probability p in [0,1], 
% rFF=1;   %r is the ratio between outlinks and inlinks selected at every "back burn" step
%k_per = k;

file_name = strcat(topology,'.csv');
adj = csvread(file_name);
n = size(adj);
n1 = n(1);
n2 = n(2);
adj_x = adj(1:n1/2, 1:n2/2);
adj_y = adj(n1/2+1:n1 , n2/2+1:n2);

adj_xy= adj(n1/2+1:n1,1:n2/2);
adj_yx= adj(1:n1/2,n2/2+1:n2);
n1 = n1/2;
n2 = n2/2;

coord=zeros(1,2);
for i=1:n1
    x=rand(1);
    y=rand(1);
    coord=[coord;x y];
end
coord=coord(2:length(coord),:);
NF_period=zeros(1,1);
%period=[2 5 10 15 20 25 40 50 75 100];
%period=[2 4 5 7 10 12 15 20 25 40 50 75 100];
%period=[5];
period=[200];   %period of time that we want to look into the propagation process

NF=zeros(1,length(period));
%for time=1:length(period)

vector_n_fail_time = cell(1,100);
vector_failed_low_betx=cell(1,100);
vector_failed_low_bety=cell(1,100);
vector_failed_low_degx=cell(1,100);
vector_failed_low_degy=cell(1,100);
vector_failed_high_betx=cell(1,100);
vector_failed_high_bety=cell(1,100);
vector_failed_high_degx=cell(1,100);
vector_failed_high_degy=cell(1,100);

vector_n_fail_time2 = cell(1,100);
vector_failed_low_betx2=cell(1,100);
vector_failed_low_bety2=cell(1,100);
vector_failed_low_degx2=cell(1,100);
vector_failed_low_degy2=cell(1,100);
vector_failed_high_betx2=cell(1,100);
vector_failed_high_bety2=cell(1,100);
vector_failed_high_degx2=cell(1,100);
vector_failed_high_degy2=cell(1,100);

vector_n_fail_time3 = cell(1,100);
vector_failed_low_betx3=cell(1,100);
vector_failed_low_bety3=cell(1,100);
vector_failed_low_degx3=cell(1,100);
vector_failed_low_degy3=cell(1,100);
vector_failed_high_betx3=cell(1,100);
vector_failed_high_bety3=cell(1,100);
vector_failed_high_degx3=cell(1,100);
vector_failed_high_degy3=cell(1,100);

vector_n_fail_time4 = cell(1,100);
vector_failed_low_betx4=cell(1,100);
vector_failed_low_bety4=cell(1,100);
vector_failed_low_degx4=cell(1,100);
vector_failed_low_degy4=cell(1,100);
vector_failed_high_betx4=cell(1,100);
vector_failed_high_bety4=cell(1,100);
vector_failed_high_degx4=cell(1,100);
vector_failed_high_degy4=cell(1,100);

vector_n_fail_time7 = cell(1,100);
vector_failed_low_betx7=cell(1,100);
vector_failed_low_bety7=cell(1,100);
vector_failed_low_degx7=cell(1,100);
vector_failed_low_degy7=cell(1,100);
vector_failed_high_betx7=cell(1,100);
vector_failed_high_bety7=cell(1,100);
vector_failed_high_degx7=cell(1,100);
vector_failed_high_degy7=cell(1,100);

vector_n_fail_time9 = cell(1,100);
vector_failed_low_betx9=cell(1,100);
vector_failed_low_bety9=cell(1,100);
vector_failed_low_degx9=cell(1,100);
vector_failed_low_degy9=cell(1,100);
vector_failed_high_betx9=cell(1,100);
vector_failed_high_bety9=cell(1,100);
vector_failed_high_degx9=cell(1,100);
vector_failed_high_degy9=cell(1,100);

vector_S_time1=cell(1,100);
vector_S_time2=cell(1,100);
vector_S_time3=cell(1,100);
vector_S_time4=cell(1,100);
vector_S_time7=cell(1,100);
vector_S_time9=cell(1,100);

vector_Dgr1_11=cell(1,100);
vector_Dgr1_12=cell(1,100);
vector_Dgr2_11=cell(1,100);
vector_Dgr2_12=cell(1,100);

vector_Dgr1_21=cell(1,100);
vector_Dgr1_22=cell(1,100);
vector_Dgr2_21=cell(1,100);
vector_Dgr2_22=cell(1,100);

vector_Dgr1_112=cell(1,100);
vector_Dgr1_122=cell(1,100);
vector_Dgr2_112=cell(1,100);
vector_Dgr2_122=cell(1,100);

vector_F1=cell(1,100);
vector_F2=cell(1,100);
vector_F3=cell(1,100);
vector_F4=cell(1,100);
vector_F7=cell(1,100);
vector_F9=cell(1,100);




for timer=1:100
    
% Erdos Renyi Model/Ramdon Network
%            adj_x = random_graph(n1,p1);   %creates the adjacancy matrix for the nodes inside network X
%            adj_y = random_graph(n2,p2);   %creates the adjacancy matrix for the nodes inside network Y
% 

% Barabasi Albert Model/Scale Free Network
%            seed =[0 1 0 0 1;1 0 0 1 0;0 0 0 1 0;0 1 1 0 0;1 0 0 0 0];  
% %            adj_x = SFNG(n1, m1, seed); 
%             adj_y = SFNG(n2, m2, seed); 

%Watts Strogatz Model/Small Workd Network
%         adj_x = small_world_graph(n1, deg1, rp1);
%          adj_y = small_world_graph(n2, deg2, rp2);







%Forest fire model
% T1=n1;
% T2=n2;
% L1 = forestFireModel(T1,pFF,rFF);
% L2 = forestFireModel(T2,pFF,rFF);
% adj_x=zeros(T1,T1);
% adj_y=zeros(T2,T2);
% for i=1:T1
%     g=L1{i};
%     for j=1:length(g)
%     adj_x(i,g(j))=1;
%     adj_x(g(j),i)=1;
%     end
% end
% 
% for i=1:T2
%     g=L2{i};
%     for j=1:length(g)
%     adj_y(i,g(j))=1;
%     adj_y(g(j),i)=1;
%     end
% end

%*****checking self loops in adj matrix
% ghotr1=zeros(1,1);
% ghotr2=zeros(1,1);
% for i=1:length(adj_x)
%     if adj_x(i,i)==1
%         ghotr1=[ghotr1 i];
%     end
%     if adj_y(i,i)==1
%         ghotr2=[ghotr2 i];
%     end
% end
%*********************
%*********************
%gplot(adj_x,coord,'-o');
degree_net_x=zeros(1,1);
for i=1:n1
    degree_node_i=sum(adj_x(:,i));
    degree_net_x=[degree_net_x degree_node_i];
end
degree_net_x=degree_net_x(1,2:length(degree_net_x));  %shows the degree of each node in network X(if a node has one parent its degree is one)
min_degree_net_x=min(degree_net_x);    %calculates the minimum degree of the nodes(minimum number of parnets of a node)
max_degree_net_x=max(degree_net_x);
index_max_degree_x=find(degree_net_x==max_degree_net_x);  
ND=unique(degree_net_x);
ND1=sort(ND,'descend');
node_degree_order_x=zeros(1,1);
for e=1:length(ND1)
    Q=find(degree_net_x==ND1(1,e));
    node_degree_order_x=[node_degree_order_x Q];
end
 node_degree_order_x=node_degree_order_x(1,2:length(node_degree_order_x)); %the index of the nodes of network X are sorted from nodes with highest degree to lowest
 
 degree_net_y=zeros(1,1);
for i=1:n2
    degree_node_i=sum(adj_y(:,i));
    degree_net_y=[degree_net_y degree_node_i];
end
degree_net_y=degree_net_y(1,2:length(degree_net_y));  %shows the degree of each node in network Y(if a node has one parent its degree is one)
dgrx_srted= sort(degree_net_x);
dgry_srted= sort(degree_net_y);
%*********************************
%********The following lines are finding the nodes with highest degree, and
%********highest centrality
%10 high percentage of the degree
wx=dgrx_srted(1,91);    %since we have 100 nodes, the upper 10% would be the last 10 degrees
wy=dgry_srted(1,91);
aqx=unique(dgrx_srted);
aqy=unique(dgry_srted);
f1=find(aqx==wx);
f2=find(aqy==wy);
T1=f1:length(aqx);
T2=f2:length(aqy);
node_index_high_10_degree=zeros(1,1);
for i=1:length(T1)
    r=T1(i);
    g=aqx(r);
    node_index_high_10_degree=[node_index_high_10_degree find(degree_net_x==g)];
end
node_index_high_10_degree=node_index_high_10_degree(1,2:length(node_index_high_10_degree)); %shows the index of the nodes in network X with 10% of highest degree

node_indey_high_10_degree=zeros(1,1);
for k=1:length(T2)
    r=T2(k);
    g=aqy(r);
    node_indey_high_10_degree=[node_indey_high_10_degree find(degree_net_y==g)];
end
node_indey_high_10_degree=node_indey_high_10_degree(1,2:length(node_indey_high_10_degree)); %shows the index of the nodes in Y that are in 10% of the highest degree
%****************************************
%10%low percentage degree
wx=dgrx_srted(1,10);
wy=dgry_srted(1,10);
aqx=unique(dgrx_srted);
aqy=unique(dgry_srted);
f1=find(aqx==wx);
f2=find(aqy==wy);
T1=1:f1;
T2=1:f2;
node_index_low_10_degree=zeros(1,1);
for i=1:length(T1)
    r=T1(i);
    g=aqx(r);
    node_index_low_10_degree=[node_index_low_10_degree find(degree_net_x==g)];
end
node_index_low_10_degree=node_index_low_10_degree(1,2:length(node_index_low_10_degree)); %shows the index of the nodes in network X with 10% of highest degree

node_indey_low_10_degree=zeros(1,1);
for j=1:length(T2)
    r=T2(j);
    g=aqy(r);
    node_indey_low_10_degree=[node_indey_low_10_degree find(degree_net_y==g)];
end
node_indey_low_10_degree=node_indey_low_10_degree(1,2:length(node_indey_low_10_degree)); %shows the index of the nodes in Y that are in 10% of the highest degree

%*****************************************
min_degree_net_y=min(degree_net_y);    %calculates the minimum degree of the nodes(minimum number of parnets of a node)
max_degree_net_y=max(degree_net_y);
index_max_degree_y=find(degree_net_y==max_degree_net_y);
index_max_degree_y=n1+index_max_degree_y;
ND0=unique(degree_net_y);
ND2=sort(ND0,'descend');
node_degree_order_y=zeros(1,1);
for e=1:length(ND2)
    Q=find(degree_net_y==ND2(1,e));
    node_degree_order_y=[node_degree_order_y Q];
end                                                                             %Samuel- X or Y ?
 node_degree_order_y=node_degree_order_y(1,2:length(node_degree_order_y)); %the index of the nodes of network X are sorted from nodes with highest degree to lowest
 
 

%******************************
initial_GC_X=length(largestcomponent(adj_x));   %Giant component of the graphX
initial_GC_Y=length(largestcomponent(adj_y));   %Giant component of the graphY





%***creating random initial failure
k=floor(k_per/100 * n1)
% F=ceil(n1*rand(1,90)); 
% F=unique(F);
F = randperm(n1);
F1=F(1,1:2*k);     % random initial failures in network A
%-----
% a = n1;
% b = n1+n2;
% r = ceil((b-a)*rand(1,90)+a );
% F=unique(r);
F= randperm(n2) +n1;
F7=F(1:2*k);    % random initial failures in network B  -added
%--------
% a = n1;
% b = n1+n2;
% r = ceil((b-a)*rand(1,90)+a );
% F=unique(r);
F= randperm(n2) +n1;
F_b=F(1:k);

% F=ceil(n1*rand(1,90));
% F=unique(F);
F0F = randperm(n1);
F_a=F(1,1:k);
F2= [ F_a F_b ];   % random initial failures in network A and B

vector_F1{timer}=F1;
vector_F2{timer}=F2;

%------------- % select nodes of A or/and B with higest centrality  
F3 = node_degree_order_x(1,1:2*k);
F9=node_degree_order_y(1,1:2*k);  %  -added
F9=F9+n1;

F_a=node_degree_order_x(1,1:k);
F_b=node_degree_order_y(1,1:k);
F_b=F_b+n1;
F4 = [ F_a F_b ];

vector_F3{timer}=F3;
vector_F4{timer}=F4;

vector_F7{timer}=F7;
vector_F9{timer}=F9;



pmax=0.8;  %the prob that a node fails if all parents are failed
kx=0.3;    %Threshold for nodes inside network X(minimum fraction of neighbors that should be failed before a node can fail)
ky=0.3;    %Threshold for nodes inside network Y
kxy=0.5;   %Threshold for parents of nodes of network Y, that are member of network X
kyx=0.5;   %Threshold for parents of nodes of network X, that are member of network Y
Ax=adj_x;
Ay=adj_y;
Axy=adj_xy;
Ayx=adj_yx;
GCX=zeros(1,1);
GCY=zeros(1,1);
 
 for scenario = 1:6
     if scenario==1
          F0=F1;
     elseif scenario==2
         F0=F2;
     elseif scenario==3
         F0=F3;
     elseif scenario==4
         F0=F4;
     elseif scenario == 5
         F0=F7;
     elseif scenario ==6
         F0=F9;
     end 
     
 for time=1:length(period)
     T=period(1,time);

[N_failed,V_state,S_time] = cascading_failure_fraction_last(F0,Ax,Ay,Axy,Ayx,kx,ky,kxy,kyx,pmax,T);

   
number_of_failed(1,time)=N_failed;
for i=1:period(1)
    f1=S_time{i}; %this shows the state of network at time step i
    f2=find(f1==1);
%     for j=1:length(F0)
%         yy=F0(j);
%         yy1=find(f2==yy);
%         f2(yy1)=[];
%     end
    list_failed{i}=f2;
end
%******We generate GML files for diff. time steps that is used to have
%******colored graphs
%generating gml files over time
% for qq=1:period(1)
%     indx=qq;
% adjacencyx=adj_x;
% adjacencyy=adj_y;
% adjacencyxy=adj_xy;
% adjacencyyx=adj_yx;
% failed=list_failed{qq};
% graphtogml(indx, adjacencyx, adjacencyy, adjacencyyx, adjacencyxy, failed)
% end
%**************************
n_fail_time=zeros(1,1);
 for k=1:period(1)
     g=length(list_failed{k});
     n_fail_time=[n_fail_time g];
 end
    n_fail_time=n_fail_time(1,2:length(n_fail_time));
    
    
    failed_low_degx=zeros(1,period(1));
    failed_low_degy=zeros(1,period(1));
    for k=1:period(1)
        a=list_failed{k};
        sum_degx=0;
        sum_degy=0;
        for i=1:length(node_index_low_10_degree)
            c1=length(find(a==node_index_low_10_degree(i)));
            sum_degx=sum_degx+c1;
        end
        for j=1:length(node_indey_low_10_degree)
            c2=length(find(a==100+node_indey_low_10_degree(j)));
            sum_degy=sum_degy+c2;
        end
        
        failed_low_degx(k)=sum_degx;
        failed_low_degy(k)=sum_degy;
    end
       failed_low_degx=(100/length(node_index_low_10_degree)).*failed_low_degx; %shows the percentage of nodes of network X with low degree that are failed over time
       failed_low_degy=(100/length(node_indey_low_10_degree)).*failed_low_degy; %shows the percentage of nodes of network Y with low degree that are failed over time
       
    
    failed_high_degx=zeros(1,period(1));
    failed_high_degy=zeros(1,period(1));
    for k=1:period(1)
        a=list_failed{k};
        sum_degx=0;
        sum_degy=0;
        for j=1:length(node_index_high_10_degree)
            c1=length(find(a==node_indexvector_n_fail_time_high_10_degree(j)));
            sum_degx=sum_degx+c1;
        end
        for z=1:length(node_indey_high_10_degree)
            c2=length(find(a==100+node_indey_high_10_degree(z)));
            sum_degy=sum_degy+c2;
        end
        
        failed_high_degx(k)=sum_degx;
        failed_high_degy(k)=sum_degy;
    end
       failed_high_degx=(100/length(node_index_high_10_degree)).*failed_high_degx; %shows the percentage of nodes of network X with high degree that are failed over time
       failed_high_degy=(100/length(node_indey_high_10_degree)).*failed_high_degy; %shows the percentage of nodes of network Y with high degree that are failed over time

%**************Here we make some plots for failed nodes with high/low centrality,etc

if scenario ==1
    vector_n_fail_time{timer} = n_fail_time; 
%     vector_failed_low_betx{timer}=failed_low_betx;
%     vector_failed_low_bety{timer}=failed_low_bety;
    vector_failed_low_degx{timer}=failed_low_degx;
     vector_failed_low_degy{timer}=failed_low_degy;
%     vector_failed_high_betx{timer}=failed_high_betx;
%     vector_failed_high_bety{timer}=failed_high_bety;
    vector_failed_high_degx{timer}=failed_high_degx;
    vector_failed_high_degy{timer}=failed_high_degy;
    
elseif  scenario ==2
    vector_n_fail_time2{timer} = n_fail_time; 
%     vector_failed_low_betx2{timer}=failed_low_betx;
%     vector_failed_low_bety2{timer}=failed_low_bety;
    vector_failed_low_degx2{timer}=failed_low_degx;
    vector_failed_low_degy2{timer}=failed_low_degy;
%     vector_failed_high_betx2{timer}=failed_high_betx;
%     vector_failed_high_bety2{timer}=failed_high_bety;
    vector_failed_high_degx2{timer}=failed_high_degx;
    vector_failed_high_degy2{timer}=failed_high_degy;

elseif scenario ==3
    vector_n_fail_time3{timer} = n_fail_time; 
%     vector_failed_low_betx3{timer}=failed_low_betx;
%     vector_failed_low_bety3{timer}=failed_low_bety;
    vector_failed_low_degx3{timer}=failed_low_degx;
    vector_failed_low_degy3{timer}=failed_low_degy;
%    vector_failed_high_betx3{timer}=failed_high_betx;
%    vector_failed_high_bety3{timer}=failed_high_bety;
    vector_failed_high_degx3{timer}=failed_high_degx;
    vector_failed_high_degy3{timer}=failed_high_degy;
    
elseif scenario==4
     vector_n_fail_time4{timer} = n_fail_time; 
%     vector_failed_low_betx4{timer}=failed_low_betx;
%     vector_failed_low_bety4{timer}=failed_low_bety;
    vector_failed_low_degx4{timer}=failed_low_degx;
    vector_failed_low_degy4{timer}=failed_low_degy;
%     vector_failed_high_betx4{timer}=failed_high_betx;
%     vector_failed_high_bety4{timer}=failed_high_bety;
    vector_failed_high_degx4{timer}=failed_high_degx;
    vector_failed_high_degy4{timer}=failed_high_degy;
elseif scenario==5
    vector_n_fail_time7{timer} = n_fail_time; 
%     vector_failed_low_betx7{timer}=failed_low_betx;
%     vector_failed_low_bety7{timer}=failed_low_bety;
    vector_failed_low_degx7{timer}=failed_low_degx;
    vector_failed_low_degy7{timer}=failed_low_degy;
%     vector_failed_high_betx7{timer}=failed_high_betx;
%     vector_failed_high_bety7{timer}=failed_high_bety;
    vector_failed_high_degx7{timer}=failed_high_degx;
    vector_failed_high_degy7{timer}=failed_high_degy;
elseif scenario==6
    vector_n_fail_time9{timer} = n_fail_time; 
%     vector_failed_low_betx9{timer}=failed_low_betx;
%     vector_failed_low_bety9{timer}=failed_low_bety;
    vector_failed_low_degx9{timer}=failed_low_degx;
    vector_failed_low_degy9{timer}=failed_low_degy;
%     vector_failed_high_betx9{timer}=failed_high_betx;
%     vector_failed_high_bety9{timer}=failed_high_bety;
    vector_failed_high_degx9{timer}=failed_high_degx;
    vector_failed_high_degy9{timer}=failed_high_degy;
    
end

failed=find(V_state==1);
failed_id{time}=failed;
loct=zeros(1,1);
for u=1:length(F0)
loct=[loct find(failed==F0(u))];
end
loct=loct(1,2:length(loct));
failed(loct)=[];  %this shows the id of failed nodes(withut the initial failure)

%***caluclation of the degree of the failed nodes
Dgr1=zeros(1,1);
Dgr2=zeros(1,1);
for u=1:length(failed)
    if failed(u)<=n1
        
        Dgr1=[Dgr1 degree_net_x(failed(u))];
        adj_x(failed(u),:)=0;
        adj_x(:,failed(u))=0;
    else
        Dgr2=[Dgr2 degree_net_y(failed(u)-n1)]; 
        adj_y(failed(u)-n1,:)=0;
        adj_y(:,failed(u)-n1)=0;
    end
end
  Dgr1=Dgr1(2:length(Dgr1)); %degree of failed nodes in network X
  Dgr2=Dgr2(2:length(Dgr2)); %degree of failed nodes in network Y
  failed_x_period{time}=Dgr1;  %degree of the failed nodes in X for period of time
%%%%%%%%%%%%%%%%%%%%%%%%%%WE are saving the results
% save initial_fail_set_ER1 initial_fail_set 
% save avr_ER1 avr
% save degree_dis_X_ER1 degree_dis_X
% save degree_dis_Y_ER1 degree_dis_Y
% 
% save failed_x_period_counter_ER1 failed_x_period_counter
% save failed_y_period_counter_ER1 failed_y_period_counter
% save fxdegree_in_5_step_ER1 fxdegree_in_5_step
% save fxdegree_in_10_step_ER1 fxdegree_in_10_step
% save fxdegree_in_20_step_ER1 fxdegree_in_20_step
% save fxdegree_in_50_step_ER1 fxdegree_in_50_step
% 
% save fydegree_in_5_step_ER1 fydegree_in_5_step
% save fydegree_in_10_step_ER1 fydegree_in_10_step
% save fydegree_in_20_step_ER1 fydegree_in_20_step
% save fydegree_in_50_step_ER1 fydegree_in_50_step
  failed_y_period{time}=Dgr2;
  GcompX=largestcomponent(adj_x);
  GcompY=largestcomponent(adj_y);
  GCX=[GCX length(GcompX)];
  GCY=[GCY length(GcompY)];
  
  if scenario==1
    vector_Dgr1_11{timer}=Dgr1;
    vector_Dgr1_12{timer}=Dgr2;
    vector_S_time1{timer}=S_time;
  
  elseif scenario==2
    vector_Dgr1_21{timer}=Dgr1;
    vector_Dgr1_22{timer}=Dgr2;
    vector_S_time2{timer}=S_time;
  elseif scenario==3
    vector_Dgr2_11{timer}=Dgr1;
    vector_Dgr2_12{timer}=Dgr2;
    vector_S_time3{timer}=S_time;
  elseif scenario==4
    vector_Dgr2_21{timer}=Dgr1;
    vector_Dgr2_22{timer}=Dgr2;
    vector_S_time4{timer}=S_time;
  elseif scenario==5
    vector_Dgr1_112{timer}=Dgr1;
    vector_Dgr1_122{timer}=Dgr2;
    vector_S_time7{timer}=S_time;
  elseif scenario==6
    vector_Dgr2_112{timer}=Dgr1;
    vector_Dgr2_122{timer}=Dgr2;
    vector_S_time9{timer}=S_time;

  end
 end
    GCX=GCX(1,2:length(GCX));  %the size of giant component of network X for diff period of time
    GCY=GCY(1,2:length(GCY));  %the size of giant component of network Y for diff period of time
    GC_X{timer}=GCX;
    GC_Y{timer}=GCY;
    initial_GCX{timer}=initial_GC_X;
    initial_GCY{timer}=initial_GC_Y;
    
NF=[NF;number_of_failed];
failed_x_period_counter{timer}=failed_x_period; %each cell shows the degree of the failed nodes for one topology and during diff periods of time
failed_y_period_counter{timer}=failed_y_period; %each cell shows the degree of the failed nodes for one topology and during diff periods of time
degree_dis_X{timer}=degree_net_x;  %shows the degree of each node for network X
degree_dis_Y{timer}=degree_net_y;  %shows the degree of each node for network Y
%NF=mean(number_of_failed);
%NF_period=[NF_period NF];
initial_fail_set{timer}=F0;

 
 end
end

%save scenari1 failed_x_period_counter failed_y_period_counter degree_dis_X degree_dis_Y initial_fail_set


disp(k*100/n1);
k_per = num2str(k_per);

temp = strcat("scenari-11&12","-",topology,"k=",k_per,"%.mat") ;
save(temp, "vector_n_fail_time" ,"vector_failed_low_degx", "vector_failed_low_degy" , "vector_failed_high_degx" ,"vector_failed_high_degy" , "vector_Dgr1_11", "vector_Dgr1_12" ,"vector_F1", "vector_F2", "vector_n_fail_time2" , "vector_failed_low_degx2", "vector_failed_low_degy2", "vector_failed_high_degx2", "vector_failed_high_degy2", "vector_Dgr1_21", "vector_Dgr1_22", "vector_S_time1", "vector_S_time2", "vector_n_fail_time7", "vector_failed_low_degy7", "vector_failed_high_degx7", "vector_failed_high_degy7", "vector_Dgr1_112", "vector_Dgr1_122", "vector_F7", "vector_S_time7")
temp = strcat("scenari-21&22","-",topology,"k=",k_per,"%.mat") ;
save( temp ,"vector_n_fail_time3" , "vector_failed_low_degx3" ,"vector_failed_low_degy3", "vector_failed_high_degx3", "vector_failed_high_degy3", "vector_Dgr2_11", "vector_Dgr2_12", "vector_F3", "vector_F4", "vector_n_fail_time4", "vector_failed_low_degx4", "vector_failed_low_degy4", "vector_failed_high_degx4", "vector_failed_high_degy4", "vector_Dgr2_21", "vector_Dgr2_22", "vector_S_time3" ,"vector_S_time4",  "vector_n_fail_time9" ,  "vector_failed_low_degx9", "vector_failed_low_degy9", "vector_failed_high_degx9", "vector_failed_high_degy9", "vector_Dgr2_112", "vector_Dgr2_122", "vector_F9", "vector_S_time9")          

 end