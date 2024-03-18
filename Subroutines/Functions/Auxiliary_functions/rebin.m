function [x_grid_new, y_grid_new] = rebin(x_grid_new, x_grid_orig, y_grid_orig)
delta = x_grid_new(2)-x_grid_new(1);
y_grid_new = zeros(size(x_grid_new));
ncounter = zeros(size(x_grid_new));

for i = 1:length(x_grid_orig)
   for j=2:length(x_grid_new)
      if x_grid_orig(i) >= x_grid_new(j-1)-delta/2 && x_grid_orig(i) < x_grid_new(j)-delta/2
          y_grid_new(j-1)  = y_grid_new(j-1) + y_grid_orig(i);
          ncounter(j-1) = ncounter(j-1) + 1;
      end
   end
end
y_grid_new = y_grid_new./ncounter;
x_grid_new(isnan(y_grid_new)) = [];
y_grid_new(isnan(y_grid_new)) = [];
end