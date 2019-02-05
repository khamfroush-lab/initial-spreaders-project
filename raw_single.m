n = 492;
clear all;
%name = 'raw_data-scenari-11&12n=100sparse&random-ER-SWk=5%std';
name = 'scenari-11&12-HIV-Human-Genetick=20%';

load(name) ;
avg_n_fail_time=0;
avg_n_fail_time2=0;



% if(size(vector_adjx{1})~=n/2)
%     msg = 'Adjx matrix not valid-';
%     s_ad = size(vector_adjx{1});
%     s_ad = num2str(s_ad);
%     disp(size(vector_adjx{1}));
%     
%     error(strcat(msg,s_ad)) ;
%     
%     
% end


for i = 1:200
    sum=0;
    
    for j =1:100
        a= vector_n_fail_time{j};
        
        sum = sum + a(i);
        
    end
    avg_n_fail_time=[avg_n_fail_time (sum/100) ];
    
end
avg_n_fail_time=avg_n_fail_time(2:length(avg_n_fail_time));


vector1_25 = 0 ;
vector1_75= 0;
vector2_25 =0;
vector2_75 = 0;


for i = 1:200
    sum=0;
    entry = zeros(1,100);
    for j =1:100
        c= vector_n_fail_time{j};
        
        entry(j) = c(i);
        
    end
    % entry = sort(entry);
    S = std(entry);
    
    vector1_25=[vector1_25 (avg_n_fail_time(i)-S)];
    vector1_75=[ vector1_75 (avg_n_fail_time(i)+S) ];
    
    
end

vector1_25 = vector1_25(2:length(vector1_25));
vector1_75 = vector1_75(2:length(vector1_75));




for i = 1:200
    
    sum2=0;
    for j =1: 100
        
        b = vector_n_fail_time2{j};
        
        sum2 = sum2 + b(i);
    end
    
    avg_n_fail_time2=[avg_n_fail_time2 (sum2/100) ];
    
end
avg_n_fail_time2=avg_n_fail_time2(2:length(avg_n_fail_time2));


for i = 1:200
    sum=0;
    entry = zeros(1,100);
    for j =1:100
        d= vector_n_fail_time2{j};
        
        entry(j) = d(i);
        
    end
    % entry = sort(entry);
    S = std(entry);
    
    vector2_25=[vector2_25 (avg_n_fail_time2(i)-S) ];
    vector2_75=[ vector2_75 (avg_n_fail_time2(i)+S)];
    
    
end

vector2_25 = vector2_25(2:length(vector2_25));
vector2_75 = vector2_75(2:length(vector2_75));


matrix_csv = [ avg_n_fail_time ; vector1_25 ;vector1_75; avg_n_fail_time2 ; vector2_25 ; vector2_75 ] ;

csvwrite(strcat(name,'.csv'),matrix_csv);

