function i_bnd = boundary2(coordinates,A)

% A = full(A);
% 
% % delaunay triangulation uses nodal coordinates to fill in complex hull
% % with triangles
% try
%     triangles = delaunay(coordinates);
% catch
%     triangles = delaunay(coordinates+...
%         randn(size(coordinates))*max(max(coordinates))*1e-12);    
% end
%     
% n_tri = size(triangles,1);
% 
% TR = triangulation(triangles,coordinates);
% i_bnd = freeBoundary(TR);
% 
% % figure(12322);clf;
% % triplot(TR,'b');hold on
% % xlocs = coordinates(:,1);
% % ylocs = coordinates(:,2);
% % for i = 1:size(i_bnd,1)
% %     plot(xlocs(i_bnd(i,:)),ylocs(i_bnd(i,:)),'o-','linewidth',2)
% % %     pause
% % end
% % % 
% % pause
% 
% 
% % i_cut_triangle = true;
% % 
% 
% 
% % check connectivity of boundary lines
% linear_ind = sub2ind(size(full(A)),i_bnd(:,1),i_bnd(:,2));
% i_delete = find(~(A(linear_ind)));
% i_delete = i_delete(1);
% while ~isempty(i_delete)
%         
%     i_bnd_add = [];
%     for i = 1:length(i_delete)
% %         i_bnd(i_delete(i),:)
%         i_cut = find(any(triangles == i_bnd(i_delete(i),1),2) & ...
%                          any(triangles == i_bnd(i_delete(i),2),2));
%                      
% %           triangles(i_cut,:)
%         if ~isempty(i_cut)
%             j_cut = find(triangles(i_cut,:) == i_bnd(i_delete(i),1) | ...
%                          triangles(i_cut,:) == i_bnd(i_delete(i),2));
%             j_not_cut = [1,2,3];
%             j_not_cut(j_cut) = [];
%         end
%         
%         i_bnd_add = [i_bnd_add;...
%                      triangles(i_cut,[j_not_cut,j_cut(1)]);...
%                      triangles(i_cut,[j_not_cut,j_cut(2)])];
%         triangles(i_cut,:) = []; 
% %         pause
%     end
%     
% %     triangles(i_cut_triangle,:) = [];
% %     n_tri = size(triangles,1);
% 
% 
% %     figure(12322);clf;
% %     triplot(TR,'b');hold on
% %     xlocs = coordinates(:,1);
% %     ylocs = coordinates(:,2);
% %     plot(xlocs(i_bnd)',ylocs(i_bnd)','ko-','linewidth',2)
% %     plot(xlocs(i_bnd(i_delete,:))',ylocs(i_bnd(i_delete,:))','rx:','linewidth',2)
% %     plot(xlocs(i_bnd_add)',ylocs(i_bnd_add)','g^--','linewidth',2)
% %     for i = 1:size(i_bnd,1)
% %         plot(xlocs(i_bnd(i,:)),ylocs(i_bnd(i,:)),'o-','linewidth',2)
% % %         pause
% %     end
% % %     
% %     pause
% 
% 
% %     TR = triangulation(triangles,coordinates);
% %     i_bnd = freeBoundary(TR);
%     i_bnd(i_delete,:) = [];
%     i_bnd = [i_bnd;i_bnd_add];
% 
% 
%     % check connectivity of boundary lines
%     linear_ind = sub2ind(size(full(A)),i_bnd(:,1),i_bnd(:,2));
%     i_delete = find(~(A(linear_ind)));
%     if ~isempty(i_delete)
%         i_delete = i_delete(1);
%     end
% end
% 
% i_bnd = unique(sort(i_bnd,2),'rows','stable');
% 
% % figure(12322);clf;
% %     triplot(TR,'b');hold on
% %     xlocs = coordinates(:,1);
% %     ylocs = coordinates(:,2);
% %     for i = 1:size(i_bnd,1)
% %         plot(xlocs(i_bnd(i,:))',ylocs(i_bnd(i,:))','o-','linewidth',2)
% %             pause
% %     end
% %     plot(xlocs(i_bnd)',ylocs(i_bnd)','o-','linewidth',2)
% %         pause
% %     end
% 
% 
% % connect segments
% i_bnd_save = i_bnd;
% i_bnd2 = i_bnd(1,:)';
% i_bnd(1,:) = nan;
% 
% i_r1 = 1;i_r2 = 1;
% while any(~isnan(i_bnd(:)))
%     
%     % see if any segments connect from beginning
%     [i_r1,i_c1] =  find(i_bnd == i_bnd2(1));
%     
%     if ~isempty(i_r1);
%         i_bnd2 = [i_bnd(i_r1(1),(rem(i_c1(1),2))+1); i_bnd2];
%         i_bnd(i_r1,:) = nan;
%     end
%     
%     % see if any segments connect from end
%     [i_r2,i_c2] =  find(i_bnd == i_bnd2(end));
%     
%     if ~isempty(i_r2)        
%         i_bnd2 = [i_bnd2;i_bnd(i_r2(1),(rem(i_c2(1),2))+1)];
%         i_bnd(i_r2(1),:) = nan;        
%     end
%     
%     if (isempty(i_r1) & isempty(i_r2))
%         i_bnd2 = i_bnd_save;
%         break
%     end
%     
% %     
% %     figure(12323);clf;
% %     triplot(TR,'b');hold on
% %     xlocs = coordinates(:,1);
% %     ylocs = coordinates(:,2);
% %     plot(xlocs(i_bnd2)',ylocs(i_bnd2)','ko-','linewidth',2)
% %         pause
%     
% end
% 
% i_bnd = i_bnd2;
% 
% % figure(12322);clf;
% % triplot(TR,'b');hold on
% % xlocs = coordinates(:,1);
% % ylocs = coordinates(:,2);
% % % for i = 1:size(i_bnd,1)
% %     plot(xlocs(i_bnd),ylocs(i_bnd),'o-','linewidth',2)
% %         pause
% % end
% 

%% Method2

A = full(A);

% delaunay triangulation uses nodal coordinates to fill in complex hull
% with triangles
triangles = delaunay(coordinates);
n_tri = size(triangles,1);

TR = triangulation(triangles,coordinates);
i_bnd = freeBoundary(TR);

% figure(12322);clf;
% triplot(TR,'b');hold on
% xlocs = coordinates(:,1);
% ylocs = coordinates(:,2);
% plot(xlocs(i_bnd),ylocs(i_bnd),'o-','linewidth',2)
% 
% pause


i_cut_triangle = true;

while any(i_cut_triangle)
    % check connectivity of boundary lines
    linear_ind = sub2ind(size(full(A)),i_bnd(:,1),i_bnd(:,2));
    i_delete = find(~(A(linear_ind)));
    i_cut_triangle = false(n_tri,1);
    for i = 1:length(i_delete)
        i_cut_triangle = i_cut_triangle | ...
                        (any(triangles == i_bnd(i_delete(i),1),2) & ...
                         any(triangles == i_bnd(i_delete(i),2),2));
    end

    triangles(i_cut_triangle,:) = [];
    n_tri = size(triangles,1);

    TR = triangulation(triangles,coordinates);
    i_bnd = freeBoundary(TR);

%     figure(12322);clf;
%     triplot(TR,'b');hold on
%     xlocs = coordinates(:,1);
%     ylocs = coordinates(:,2);
%     plot(xlocs(i_bnd),ylocs(i_bnd),'o-','linewidth',2)
%     
%     pause(0.1)
end




% connect segments
i_bnd2 = i_bnd(1,:)';
i_bnd(1,:) = nan;
while any(~isnan(i_bnd(:)))
    
    % see if any segments connect from beginning
    [i_r1,i_c1] =  find(i_bnd == i_bnd2(1));
    
    if ~isempty(i_r1);
        i_bnd2 = [i_bnd(i_r1(1),(rem(i_c1(1),2))+1); i_bnd2];
        i_bnd(i_r1,:) = nan;
    end
    
    % see if any segments connect from end
    [i_r2,i_c2] =  find(i_bnd == i_bnd2(end));
    
    if ~isempty(i_r2)        
        i_bnd2 = [i_bnd2;i_bnd(i_r2(1),(rem(i_c2(1),2))+1)];
        i_bnd(i_r2(1),:) = nan;        
    end
end

i_bnd = i_bnd2;


