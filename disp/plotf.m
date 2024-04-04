function fig = plotf(x,y,line_feature)

if nargin == 1
    y = x;
    x = 1:length(y);
    plot(x,y,'k-','LineWidth',3); grid on;
end

if nargin == 2 
    if isnumeric(y) == 0
        line_feature = y;
        y = x;
        x = 1:length(y);
        plot(x,y,line_feature,'LineWidth',3); grid on;
    else 
        plot(x,y,'k-','LineWidth',3); grid on;
    end
end

if nargin == 3
    plot(x,y,line_feature,'LineWidth',3);grid on;
end

set(gca,'FontSize',20,'FontWeight','bold');
axis square

end