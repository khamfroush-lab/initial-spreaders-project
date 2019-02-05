n = 500;

avg_n_fail_time=0;
avg_n_fail_time2=0;

for i = 1:200
    sum2=0;
    
    for j =1:100
        a= vector_S_time1{j};
        a=a{i};
        a=a(1:n);
        b = sum(a);
        sum2 = sum2 + b;
    end
    avg_n_fail_time=[avg_n_fail_time (sum2/100) ];
end



avg_n_fail_time=avg_n_fail_time(2:length(avg_n_fail_time));

t = 1:length(avg_n_fail_time);


for i = 1:200
    sum2=0;
    for j =1: 100
        a= vector_S_time1{j};
        a=a{i};
        a=a(n+1:2*n);
        b = sum(a);
        sum2 = sum2 + b;
    end
    avg_n_fail_time2=[avg_n_fail_time2 (sum2/100) ];
end
avg_n_fail_time2=avg_n_fail_time2(2:length(avg_n_fail_time2));

figure
hold on 
plot(t,avg_n_fail_time);
xlabel('time')
ylabel('avg No. of failed nodes')

plot(t,avg_n_fail_time2);

legend(' first network ', '  second');

hold off

