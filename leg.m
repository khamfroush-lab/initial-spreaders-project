x = 1:10; y1 = 4*x;
y2 = 3*x + 5; 
figure(1), 
plot(x,y1,'b',x,y2,'r')
axP = get(gca,'Position'); 
legend('Line 1','Line 2','Location','NorthEastOutside') 
set(gca, 'Position', axP)