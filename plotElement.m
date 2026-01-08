function plotElement(elements, nodes)
    figure()
    hold on
    scatter(nodes(:,1),nodes(:,2), 'ok');
    t = 0.8;
    for ii=1:size(elements,1)
        centroid = mean(nodes(elements(ii,2:end),:),1);
        text(centroid(1),centroid(2),num2str(ii),'Color','r');
        for jj=1:(size(elements,2)-1)
            node_pos = nodes(elements(ii,jj),:);
            text_pos = t*node_pos+(1-t)*centroid;
            text(text_pos(1),text_pos(2),num2str(elements(ii,jj)),'Color','b');
        end
        pos_array = nodes(elements(ii,1:end),:);
        pos_array(end+1,:) = pos_array(1,:); %close loop
        if size(pos_array,1)>6
            temp = pos_array;
            pos_array(1,:)=temp(3,:);
            pos_array(2,:)=temp(1,:);
            pos_array(3,:)=temp(4,:);
            pos_array(4,:)=temp(2,:);
        end
        plot(pos_array(:,1),pos_array(:,2),'-k')
        
    end
    hold off
end