

%name = 'scenari-21&22-Arabidopsis-Multiplex-Genetick=10%';

%function average_1(n)
num =2;


files = dir;

disp(length(files))


for k = 3: 11
num = num +2;
load(files(k).name)
disp(files(k).name);


avg_n_fail_time=0;
avg_n_fail_time2=0;

for i = 1:200
    sum=0;
    for j =1:100
        a= vector_n_fail_time{j};
        sum = sum + a(i);
    end
    avg_n_fail_time=[avg_n_fail_time (sum/100) ];
end



avg_n_fail_time=avg_n_fail_time(2:length(avg_n_fail_time));

t = 1:length(avg_n_fail_time);


% for i = 1:200
%     sum=0;
%     for j =1: 100
%         a= vector_n_fail_time2{j};
%         sum = sum + a(i);
%     end
%     
% %     ks=sum;
% %     if( i==1)
% %         disp(ks/100);
% %     end
%     
%     avg_n_fail_time2=[avg_n_fail_time2 (sum/100) ];
% end
% avg_n_fail_time2=avg_n_fail_time2(2:length(avg_n_fail_time2));

%fig = figure('Name',name,'NumberTitle','off');

hold on 
plot(t,avg_n_fail_time);
xlabel('time')
num2 = num2str(num)
ylabel('average number of failed nodes')
%plot(t,avg_n_fail_time2);

axis([0 200 0 200]);

% fig2 = figure;
% legendflex(h2, {'Thing 1', 'Thing 2'}, 'ref', fig2)




 legend(strcat(num2,' averagedegree'));
 


hold off

%saveas(fig,strcat(name,'.pdf'))


%end


end