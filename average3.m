
files = dir;
h=0;
disp(length(files))
figure 
for k =3:2:14
   
  fprintf('%s \n', files(k).name) ;
  fprintf('%s \n', files(k+1).name)
  

  load(files(k).name);
  load(files(k+1).name);

  
  
curr_min=200;
curr_min2=200;
for i = 1: 100
    a =vector_n_fail_time3{i};
    b = vector_n_fail_time2{i};
     if min(find(a==200)) < curr_min
         curr_min= min(find(a==200));
     end
     if min(find(b==200)) < curr_min2
         curr_min2= min(find(b==200));
     end
     
     
end

limit = curr_min;
limit2=curr_min2;

avg_n_fail_time=0;
avg_n_fail_time2=0;

for i = 1:200
    sum=0;
    
    for j =1: 100
        a= vector_n_fail_time3{j};
       
        sum = sum + a(i);
       
    end
    avg_n_fail_time=[avg_n_fail_time (sum/100) ];
   

end

for i = 1:200
    
    sum2=0;
    for j =1: 100
        
        b = vector_n_fail_time2{j};
        
        sum2 = sum2 + b(i);
    end
    
    avg_n_fail_time2=[avg_n_fail_time2 (sum2/100) ];

end



avg_n_fail_time=avg_n_fail_time(2:length(avg_n_fail_time));
avg_n_fail_time2=avg_n_fail_time2(2:length(avg_n_fail_time2));



t = 1:length(avg_n_fail_time);
t2=1:length(avg_n_fail_time2);

h=h+1;
disp(h);
subplot(3,3,h);


plot(t,avg_n_fail_time);
y=ylabel('avg No. of failed nodes');
set(y,'FontSize',8);

if h==1
    %legend('SF-SF' ); 
    x=xlabel('time(SF-SF)');
    set(x,'FontSize',8);
    
elseif h==2
        %legend('ER-ER');
        x=xlabel('time(ER-ER)');
        set(x,'FontSize',8);
        
elseif h==3
       % legend('SW-SW');
        x=xlabel('time(SW-SW)');
        set(x,'FontSize',8);
        
elseif h==4
        %legend('ER-SF');
        x=xlabel('time(ER-SF)');
        set(x,'FontSize',8);
elseif h==5
    %legend('SF-SW');
    x=xlabel('time(SF-SW)');
    set(x,'FontSize',8);
elseif h==6
    %legend('ER-SW');
    x=xlabel('time(ER-SW)');
    set(x,'FontSize',8);
end

%legend(' graph1')

hold on

plot(t2,avg_n_fail_time2);

if h==3
    legend(' 6 spreaders in C' , ' 3 spreaders in each of C & P');
end

hold off


end

h=6;
%%%% and------------------

if length(files) == 14

for k =9:2:13
    
  fprintf('%s \n', files(k).name) ;
  fprintf('%s \n', files(k+1).name)
  load(files(k).name);
  load(files(k+1).name);

  disp(33333)
  
curr_min=200;
curr_min2=200;
for i = 1: 100
    a =vector_n_fail_time9{i};
    b = vector_n_fail_time2{i};
     if min(find(a==200)) < curr_min
         curr_min= min(find(a==200));
     end
     if min(find(b==200)) < curr_min2
         curr_min2= min(find(b==200));
     end
     
     
end

limit = curr_min;
limit2=curr_min2;

avg_n_fail_time=0;
avg_n_fail_time2=0;

for i = 1:200
    sum=0;
    
    for j =1: 100
        a= vector_n_fail_time9{j};
       
        sum = sum + a(i);
       
    end
    avg_n_fail_time=[avg_n_fail_time (sum/100) ];
   

end

for i = 1:200
    
    sum2=0;
    for j =1: 100
        
        b = vector_n_fail_time2{j};
        
        sum2 = sum2 + b(i);
    end
    
    avg_n_fail_time2=[avg_n_fail_time2 (sum2/100) ];

end



avg_n_fail_time=avg_n_fail_time(2:length(avg_n_fail_time));
avg_n_fail_time2=avg_n_fail_time2(2:length(avg_n_fail_time2));



t = 1:length(avg_n_fail_time);
t2=1:length(avg_n_fail_time2);

h=h+1;
disp(h);
subplot(3,3,h);


plot(t,avg_n_fail_time);
y=ylabel('avg No. of failed nodes');
set(y,'FontSize',8);

if h==7
    %legend('SF-SF' ); 
    x=xlabel('time(SF-ER)');
    set(x,'FontSize',8);
    %legend(' 6 spreaders in A' , ' 3 in A & B');
elseif h==8
        %legend('ER-ER');
        x=xlabel('time(SW-SF)');
        set(x,'FontSize',8);
elseif h==9
       % legend('SW-SW');
        x=xlabel('time(SW-ER)');
        set(x,'FontSize',8);
end

%legend(' graph1')

hold on

plot(t2,avg_n_fail_time2);



hold off


end

else
    
    
 for k =15:2:19
    
  fprintf('%s \n', files(k).name) ;
  fprintf('%s \n', files(k+1).name)
  load(files(k).name);
  load(files(k+1).name);

  disp(444)
  
curr_min=200;
curr_min2=200;
for i = 1: 100
    a =vector_n_fail_time9{i};
    b = vector_n_fail_time2{i};
     if min(find(a==200)) < curr_min
         curr_min= min(find(a==200));
     end
     if min(find(b==200)) < curr_min2
         curr_min2= min(find(b==200));
     end
     
     
end

limit = curr_min;
limit2=curr_min2;

avg_n_fail_time=0;
avg_n_fail_time2=0;

for i = 1:200
    sum=0;
    
    for j =1: 100
        a= vector_n_fail_time3{j};
       
        sum = sum + a(i);
       
    end
    avg_n_fail_time=[avg_n_fail_time (sum/100) ];
   

end

for i = 1:200
    
    sum2=0;
    for j =1: 100
        
        b = vector_n_fail_time2{j};
        
        sum2 = sum2 + b(i);
    end
    
    avg_n_fail_time2=[avg_n_fail_time2 (sum2/100) ];

end



avg_n_fail_time=avg_n_fail_time(2:length(avg_n_fail_time));
avg_n_fail_time2=avg_n_fail_time2(2:length(avg_n_fail_time2));



t = 1:length(avg_n_fail_time);
t2=1:length(avg_n_fail_time2);

h=h+1;
disp(h);
subplot(3,3,h);


plot(t,avg_n_fail_time);
y=ylabel('avg No. of failed nodes');
set(y,'FontSize',8);

if h==7
    %legend('SF-SF' ); 
    x=xlabel('time(SF-ER)');
    set(x,'FontSize',8);
    %legend(' 6 spreaders in A' , ' 3 in A & B');
elseif h==8
        %legend('ER-ER');
        x=xlabel('time(SW-SF)');
        set(x,'FontSize',8);
elseif h==9
       % legend('SW-SW');
        x=xlabel('time(SW-ER)');
        set(x,'FontSize',8);
end

%legend(' graph1')

hold on

plot(t2,avg_n_fail_time2);



hold off


end
    
    
    
end


% 

% %save image 
% 
% 
% 
% g =6;
% 
% % scenario 2 dense & random SW to SF
% 
%     
%     
%     
%     
%     
%     