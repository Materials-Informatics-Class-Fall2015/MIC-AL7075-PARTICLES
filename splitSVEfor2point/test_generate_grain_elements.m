
fid = fopen(['trial_elem_grains_' num2str(0) '.txt'],'w');

x = 10;
y = 10;
z = 10;
fprintf(fid,'x,y,z\n');
fprintf(fid, '%d,%d,%d\n', x,y,z);
for i = 1:z
    for j = 1:y
        for k = 1:x
            if(i > 5)
                fprintf(fid, '%d\n', 1);
            else
                fprintf(fid, '%d\n', 2);
            end
        end
    end
end
fclose(fid);