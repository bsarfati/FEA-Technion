function plotNodes(V)
    %Written by Ben 
    
    %Plot settings 
    bufferSize = 1;
    fontSize = 30;
    axFontSize = 30;
    markerSize = 24;
    textSpacing = 0.1;

    %Plot
    figure;
    plot(V(:,1),V(:,2), '.k','MarkerSize',markerSize);
    set(gcf,'color','w')
    for i = 1:size(V,1)
        text(V(i,1)+textSpacing,V(i,2)+textSpacing,num2str(i),'Color','b','FontSize',fontSize);
    end
    xlabel('$x$','FontSize',fontSize,'Interpreter','latex')
    ylabel('$y$','FontSize',fontSize,'Interpreter','latex')
    set(gca,'FontSize',axFontSize)

    %Set appropriate bounds to the plot
    xMax = max(V(:,1));
    yMax = max(V(:,2));
    xMin = min(V(:,1));
    yMin = min(V(:,2));
    xlim([xMin-bufferSize xMax+bufferSize])
    ylim([yMin-bufferSize yMax+bufferSize])