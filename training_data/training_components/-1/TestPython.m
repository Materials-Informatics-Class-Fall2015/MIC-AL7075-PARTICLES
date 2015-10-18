x = 20;
y = 23;
z = 20;
res = 1;
fid = fopen('test.inp','wt');
type = 'rough_crack';
options.crack_height = '1';
mesh.featureSize = '45';

% assign node sets for crack tip top and bototm
if ~isempty(strfind(type, 'crack'))
    featureSize = str2double(mesh.featureSize)*x/100;
    height = str2double(options.crack_height);
    
    if ~isempty(strfind(type, 'rough'))
        elemsRemove = ceil(str2double(options.crack_height)/100*y/res)
        center = y/2+(1-mod(y/res,2)/2-mod(elemsRemove,2)/2)*res
        crack_top = center + elemsRemove/2*res
        crack_bottom = center - elemsRemove/2*res
        crack_in = (floor(featureSize/res-.5)+1)*res
    else
        center = y/2;
        crack_top = center + height/2;
        crack_bottom = center - height/2;
        crack_in = featureSize;
    end
end

