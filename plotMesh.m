function plotMesh(E,V)
    %Written by Ben 
    
    %Plot settings 
    bufferSize = 1;
    fontSize = 30;
    axFontSize = 30;
    markerSize = 24;
    textSpacing = 0.1;

    %Plot nodes
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

    %Add element numbering and edges to plot
    hold on;
    for i = 1:length(E)
        %Add element numbering according to order in E
        centroid = mean(V(E{i},:));
        text(centroid(1),centroid(2),num2str(i),'Color','r','FontSize',fontSize);

        % %Add edges according to type of element 
        % if length(E{i})>6 %element is biquad
        %     edges = V(E{i}([2 3 4 1]),:)-V(E{i}(1:4),:)
        % else %element is quad tri
        %     plot()
        % end
            
    end

