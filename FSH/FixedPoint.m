E2_values = linspace(20, 200, 10);
Ih_values = linspace(1, 10, 10);

fixedPoints = zeros(length(E2_values), length(Ih_values));

for i = 1:length(E2_values)
    for j = 1:length(Ih_values)
        E2_input = @(t) E2_values(i);
        Ih_input = @(t) Ih_values(j);
        fixedPoints(i, j) = find_fixed_points_FSH(E2_input, Ih_input);
    end
end
