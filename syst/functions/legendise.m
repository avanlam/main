function legend_cell = legendise(factors)
%LEGENDISE : Translate the list of factor names to a LATEX list of factor
% names for the graphs

legend_cell = factors.name;
for i = 1:numel(legend_cell)
    legend_cell{i} = strcat('$$', legend_cell{i},'$$');
end
legend_cell{end} = "$$\Omega$$";

end

