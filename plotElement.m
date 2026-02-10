function plotElement(elements, nodes)
        %Implement functionality for quadratic 
        if size(elements,2) > 3
            elements = elements(:,1:3);
        end

        scatter(nodes(:,1),nodes(:,2), 'ok','MarkerFaceColor','k');
        t = 0.8;
        for ii=1:size(elements,1)
            centroid = mean(nodes(elements(ii,:),:),1);
            text(centroid(1),centroid(2),num2str(ii),'Color','r');
            % for jj=1:size(elements,2)
            %     node_pos = nodes(elements(ii,jj),:);
            %     text_pos = t*node_pos+(1-t)*centroid;
            %     text(text_pos(1),text_pos(2),num2str(elements(ii,jj)),'Color','b');
            % end
            pos_array = nodes(elements(ii,:),:);
            pos_array(end+1,:) = pos_array(1,:); %close loop
            plot(pos_array(:,1),pos_array(:,2),'-k')
            
        end
        
end