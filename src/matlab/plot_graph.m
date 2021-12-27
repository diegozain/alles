function plot_graph(nodes,c,thic)
% diego domenzain
% dec 2021
% ------------------------------------------------------------------------------
% c is Adjacency matrix
c_ = triu(c);
% nodes is a list of size nnodes × (x,y)
nnodes = size(nodes,1);
% ------------------------------------------------------------------------------
figure;
hold on;
% plots all edges
for in=1:nnodes
    in_ = find(c_(in,:));
    for inn=1:numel(in_)
        % edges
        a = [nodes(in,1) , nodes(in_(inn),1)];
        b = [nodes(in,2) , nodes(in_(inn),2)];
        plot(a,b,'k-')
        % % edges weights
        % str_ = full(c_(in,in_(inn)));
        % a = sum(a)*0.5;
        % b = sum(b)*0.5;
        % text(a,b,num2str(str_),'fontsize',20);
    end
end
% plot all nodes
plot(nodes(:,1),nodes(:,2),'k.','markersize',thic)
hold off
% axis tight
ylim([min(nodes(:,2))-0.2,max(nodes(:,2))+0.2])
xlim([min(nodes(:,1))-0.2,max(nodes(:,1))+0.2])
axis square
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Graph')
simple_figure()
end
