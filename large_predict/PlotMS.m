s = size(S_star);
%s = [30,30,30];
h=figure;
set(h,'visible','off');

hold off
for i = 1:s(1)
   for j = 1:s(2)
        for k = 1:s(3)
            if(S_star(i,j,k)==0)
               plot3(i,j,k)
               hold on
            end
        end
   end
end
set(h,'visible','on');