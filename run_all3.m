 
top = [ "ER-ER" "ER-SF" "ER-SW" "SF-SF" "SF-SW" "SF-ER" "SW-SW" "SW-ER" "SW-SF" ];
inter = [ "sparse&random" "dense&random" "sparse&designed-max_max" "sparse&designed-max_min" "sparse&designed-min_min" "dense&designed-max_max" "dense&designed-max_min" "dense&designed-min_min" ];
for i = 8:8
   for j = 6:6 
       for avg = 4:2:20
           
          failuremodel_mod3(500,500,top(i),inter(j),avg);
       end
       
   end
  
end
