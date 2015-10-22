%% Visualize difference between actual and reconstructed matrix
% From http://stackoverflow.com/questions/3942892/how-do-i-visualize-a-...
% matrix-with-colors-and-values-displayed

function visualizeDelta(mat)
m = length(mat);
imagesc(mat);            %# Create a colored plot of the matrix values
colormap(flipud(gray));  %# Change the colormap to gray (so high values are
%#   black and lower values are white)

textStrings = num2str(mat(:),'%0.2f');  %# Create string from matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
[x,y] = meshgrid(1:m);   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
    'HorizontalAlignment','center');
midValue = mean(get(gca,'CLim'));  %# Get the middle value of color range
textColors = repmat(mat(:) > midValue,1,3);  %# Choose white or black for

%#   text color of the strings so
%#   they can be easily seen over
%#   the background color
set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors

set(gca,'XTick',2:m,...                         %# Change axes tick marks
    'YTick',2:m,...
    'TickLength',[0 0]);
end