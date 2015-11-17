S_star = ones(200,200,200);
blah = numel(S_star);
seeds = 0.02*blah;
side = size(S_star);
side = side(1);
for i = 1:seeds
    ind = uint32(rand(3,1)*(side-1)+1);
    S_star(ind(1),ind(2),ind(3)) = 0;
end