
% files = dir;
% 
% for k = 3:length(files)
%     
% 
%     
% load(files(k).name);
% name = files(k).name;
 name = 'scenari-11&12n=500dense&random-ER-SWk=10%-avgdeg=8.mat';
%function average_1(n)

avg_n_fail_time=0;
avg_n_fail_time2=0;

for i = 1:200
    sum=0;
    for j =1:100
        a= vector_n_fail_time9{j};
        sum = sum + a(i);
    end
    avg_n_fail_time=[avg_n_fail_time (sum/100) ];
end



avg_n_fail_time=avg_n_fail_time(2:length(avg_n_fail_time));

t = 1:length(avg_n_fail_time);

for i = 1:200
    sum=0;
    for j =1: 100
        a= vector_n_fail_time4{j};
        sum = sum + a(i);
    end
    
%     ks=sum;
%     if( i==1)
%         disp(ks/100);
%     end
    
    avg_n_fail_time2=[avg_n_fail_time2 (sum/100) ];
end
avg_n_fail_time2=avg_n_fail_time2(2:length(avg_n_fail_time2));

fig = figure('Name',name,'NumberTitle','off');

hold on 
plot(t,avg_n_fail_time);
xlabel('time')
ylabel('avg No. of failed nodes')

h2 = plot(t,avg_n_fail_time2);
%axis([0 200 0 1000]);

% fig2 = figure;
% legendflex(h2, {'Thing 1', 'Thing 2'}, 'ref', fig2)




 lgnd = legend('|F_0| spreaders in A', '|F_0|/2 in each of A & B');
 set(lgnd,'color','none');



hold off
% 
saveas(fig,strcat(name,'.pdf'))


%end

% end
