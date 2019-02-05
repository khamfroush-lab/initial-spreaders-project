 
top = [ "ER-ER" "ER-SF" "ER-SW" "SF-SF" "SF-SW" "SF-ER" "SW-SW" "SW-ER" "SW-SF" ];
inter = [ "sparse&random" "dense&random" "sparse&designed-max_max" "sparse&designed-max_min" "sparse&designed-min_min" "dense&designed-max_max" "dense&designed-max_min" "dense&designed-min_min" ];
for i = 1:1
   for j = 6:6
       failuremodel_mod5(500,500,top(i),inter(j))
  
   end

end