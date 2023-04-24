function h = plot_area(x,y,PropValues)

% Define the x and y coordinates of the area to be filled
x_fill = [x(1) x x(end)];
y_fill = [0 y 0];

% Create an area plot with the specified vertices and color
area(x_fill, y_fill, PropValues{:});