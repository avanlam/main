function legend_cell = legendise(factors)
%LEGENDISE Summary of this function goes here
%   Detailed explanation goes here
legend_cell = factors.name;
for i = 1:numel(legend_cell)
    legend_cell{i} = strcat('$$', legend_cell{i},'$$');
end
legend_cell{end} = "$$\Omega$$";

end

