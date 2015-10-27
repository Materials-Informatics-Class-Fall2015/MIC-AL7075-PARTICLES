function [grain,ElemGrainNo,V_grn,RelLength]=Circular_Packing_CG(d_gr,sigma_fit,mu_fit,geom,n_grains, el_centroid, max_Iter,n_El,V_el,Big_grain1,Big_grain2,Big_grain3, options)
% return values are
% grain: vector of grain structs fields = x,y,z, vol, r, elem
% ElemGrainNo: vector of each element mapping to its grain number
% V_grn: vector of each grain volume
% RelLength: vector of x,y,z distances for element from its grain center

disp(' ')
disp('Begin circular packing')

% Bounds on grain placement zone
aa=0.000;
xmin=0+aa;
xmax=geom.x-aa;
ymin=0+aa;
ymax=geom.y-aa;
zmin=0+aa;
zmax=geom.z-aa;

% Grain size distribution pdf parameters
%sigma_fit=0.2;
%mu_fit=-0.1;

% Create list of spherical radii values and sort in descending order

norm_grn_vol=lognrnd(mu_fit,sigma_fit,n_grains,1); % returns an array of random numbers generated from the lognormal distribution with parameters mu and sigma, where scalars n_grains and 1 are the row and column dimensions

% Reduce volume values by factor to account for imperfect packing and regulirize by maxi value
norm_grn_vol=0.55*norm_grn_vol/(max(norm_grn_vol));

% Equivalent radii
%d_gr=8;  %diameter of coarse grain
grn_size_dist=0.5*d_gr*(norm_grn_vol).^(1/3);

% Sort radii 
grn_size_dist=sort(grn_size_dist,'descend');

cent_x=el_centroid(:,1);  %el_centroid only contains elements in CP region
cent_y=el_centroid(:,2);
cent_z=el_centroid(:,3);
ElemGrainNo=int32(zeros(n_El,1));

Used_elem=[];

%% Initial placement of spheres

for ii=1:n_grains


     if ii==1 && Big_grain1 ~= 0
    %Current grain radius
    grain(ii).r=0.5*0.8193*Big_grain1; %(divide by half to calculate radius. Multiply by 0.55^(1/3) for imperfect packing)
    %Current grain target volume
    grain(ii).vol_tar=4/3*pi*grain(ii).r^3;
    elseif ii==2 && Big_grain2 ~= 0
    %Current grain radius
    grain(ii).r=0.5*0.8193*Big_grain2; %(divide by half to calculate radius. Multiply by 0.55^(1/3) for imperfect packing)
    %Current grain target volume
    grain(ii).vol_tar=4/3*pi*grain(ii).r^3;
    elseif ii==3 && Big_grain3 ~= 0
    %Current grain radius
    grain(ii).r=0.5*0.8193*Big_grain3; %(divide by half to calculate radius. Multiply by 0.55^(1/3) for imperfect packing)
    %Current grain target volume
    grain(ii).vol_tar=4/3*pi*grain(ii).r^3;
    else
    %Current grain radius
    grain(ii).r=grn_size_dist(ii);
    %Current grain target volume
    grain(ii).vol_tar=4/3*pi*grain(ii).r^3;
    end


    count=1;

    while (count<=max_Iter)
        dup=0;


     if ii==1 && Big_grain1 ~= 0
        grain(ii).x=(xmax-xmin)/2;
        grain(ii).y=(ymax-ymin)/2;
        grain(ii).z=(zmax-zmin)/2;


     elseif ii==2 && Big_grain2 ~= 0
        grain(ii).x=(xmax-xmin)/2;
        grain(ii).y=((ymax-ymin)/2-grain(1).r)/2;
        grain(ii).z=(zmax-zmin)/2;


     elseif ii==3 && Big_grain3 ~= 0
        grain(ii).x=((xmax-xmin)/2-grain(1).r)/2 ;
        grain(ii).y=(ymax-ymin)/2;
        grain(ii).z=(zmax-zmin)/2;

     else
        % Pick random location for grain sphere
        grain(ii).x=roundn(xmax-rand()*(xmax-xmin),-3);
        grain(ii).y=roundn(ymax-rand()*(ymax-ymin),-3);
        grain(ii).z=roundn(zmax-rand()*(zmax-zmin),-3);
     end 

        %Distance from grain centroid to each element centroid
        dist=((cent_x-grain(ii).x).^2+(cent_y-grain(ii).y).^2+(cent_z-grain(ii).z).^2).^(1/2);
            
        %Find elements within grain radius 
        a= find(dist<=grain(ii).r);
                
            %Check that element set is unique (elements have not been
            %used on a previous grain)
            for kk=1:numel(a)
                b=find(a(kk)==Used_elem);
                    if numel(b)>0
                        dup=2;
                        break
                    else
                        dup=1;
                    end
                        
            end
        

        if dup==1; %Grain sphere is ok at current location
            Used_elem=[Used_elem;a]; % Add used elements to "Used" array
            for kk=1:numel(a)
                ElemGrainNo(a)=ii;                
            end
            grain(ii).vol=sum(V_el(a));
            grain(ii).elem=a;
            break
        end
        
        count=count+1;
        
    end

end

x=zeros(n_grains,1);
y=zeros(n_grains,1);
z=zeros(n_grains,1);
r=zeros(n_grains,1);
V_grn=zeros(n_grains,1);

for jj=1:n_grains
    x(jj,1)=grain(jj).x;
    y(jj,1)=grain(jj).y;
    z(jj,1)=grain(jj).z;
    r(jj,1)=grain(jj).r;
    if grain(jj).vol>0
        V_grn(jj,1)=grain(jj).vol;
    end
end

% Print out txt file of grain center cordinates:
    fid = fopen(['Grain_Centers_' num2str(options.MS_number) '.txt'], 'w'); 
for jj=1:n_grains
    fprintf(fid,[num2str(x(jj,1)) ' ' num2str(y(jj,1)) ' ' num2str(z(jj,1)) '\n']); 
end
    fclose(fid);

%% Add Elements to grains w/o elements:

% unassigned elements
no_elem=find(ElemGrainNo==0);
n_no_elem=numel(no_elem);

% grains with zero volume
no_grn=find(V_grn==0);
n_no_grn=numel(no_grn);

if n_no_grn>=1
    for ii=1:n_no_grn
        dist=((cent_x(no_elem)-x(no_grn(ii))).^2+(cent_y(no_elem)-y(no_grn(ii))).^2+(cent_z(no_elem)-z(no_grn(ii))).^2).^(1/2);
        a=find(dist==min(dist));
        temp_elem=no_elem(a);
        temp_elem=temp_elem(1);
        temp_grn=no_grn(ii);
        ElemGrainNo(temp_elem)=temp_grn;
        grain(temp_grn).elem=temp_elem;
        grain(temp_grn).vol=V_el(temp_elem,1);
        V_grn(temp_grn,1)=grain(temp_grn).vol;
        % Eliminate element from list
        Used_elem=[Used_elem;temp_elem];
        no_elem(a)=[];
    end
else
end

n_no_elem=n_no_elem-n_no_grn;


%% Grow spherical grains to non-"Used_elem"

for ii=1:n_El

    if ElemGrainNo(ii,1)>0
        % Element already used
    else
        %Distance from element centroid to edge of circular outside
        dist=((cent_x(ii,1)-x).^2+(cent_y(ii,1)-y).^2+(cent_z(ii,1)-z).^2).^(1/2)-r;
        a=find(dist==min(dist));
        b=numel(a);
	[C,a]=min(dist);
        a=a(1);
        ElemGrainNo(ii)=a;
        grain(a).elem=[grain(a).elem;ii];
        grain(a).vol=V_grn(a,1)+V_el(ii,1);
        V_grn(a,1)=grain(a).vol;
    end
end


RelLength=zeros(n_El,3);

for ii=1:n_El   
%          %Vector from grain centroid to each element centroid
           RelLength(ii,1)=(cent_x(ii)- grain(ElemGrainNo(ii,1)).x);
           RelLength(ii,2)=(cent_y(ii)-grain(ElemGrainNo(ii,1)).y);
           RelLength(ii,3)=(cent_z(ii)-grain(ElemGrainNo(ii,1)).z); 
end



disp(' ')
disp('Finish circular packing')

return
end
