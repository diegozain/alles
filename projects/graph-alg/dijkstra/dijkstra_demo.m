clear
close all
clc
% ------------------------------------------------------------------------------
addpath('src');
% ------------------------------------------------------------------------------
%     dijkstra on a graph
% ------------------------------------------------------------------------------
% number of nodes in graph
nnodes= 10;
% choose source position in
% index notation
s = 1;
% want path from s to r?
% use magic vector P. (r is in index notation)
r = fix(0.5*nnodes);
% ------------------------------------------------------------------------------
% build velocity profile
% ------------------------------------------------------------------------------
% % custom
% c = [0 7 9 0 0 14; 7 0 10 15 0 0;...
%     9 10 0 11 0 2; 0 15 11 0 6 0;...
%     0 0 0 6 0 9; 14 0 2 0 9 0];
% nnodes = 6;

% random
c = randi(5,nnodes);
c = triu(c);
% this matrix is too full,
% lets make it sparser.
% also, lets get rid of entries in the diagonal
for i_=0:nnodes-1
a = randi(nnodes,[fix(nnodes*0.7),1]);
c(i_+1,a) = 0;
c(1+i_*(nnodes+1))=0;
end
c = c + c.';
% ------------------------------------------------------------------------------
figure;
subplot(121)
fancy_imagesc(c)
colormap(rainbow2(1))
colorbar('off')
title('Adjacency matrix')
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure()
% ------------------------------------------------------------------------------
c= sparse(c);
% -----------------------------------------------------------------------------
% compute field u on all nodes and
% magic vector P
[u,P] = dijkstra(c,s);

p = ray(s,r,P);
p = full(p);
% ------------------------------------------------------------------------------
% see
% ------------------------------------------------------------------------------
thet = linspace(0,2*pi,nnodes+1);
thet = thet.';
thet(1) = [];
nodes = [cos(thet) , sin(thet)];
% ------------------------------------------------------------------------------
c_ = triu(c);

% figure;
subplot(122)
hold on
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
% plot path edges
for in=1:numel(p)-1
    % edges
    a = [nodes(p(in),1) , nodes(p(in+1),1)];
    b = [nodes(p(in),2) , nodes(p(in+1),2)];
    plot(a,b,'b-','linewidth',5)
end
% plot all nodes
plot(nodes(:,1),nodes(:,2),'k.','markersize',50)
% plot path nodes
plot(nodes(p(1),1),nodes(p(1),2),'rp','markersize',10,'markerfacecolor','r')
for ip = 2:numel(p)-1
    plot(nodes(p(ip),1),nodes(p(ip),2),'.','markersize',30)
end
plot(nodes(p(end),1),nodes(p(end),2),'wp','markersize',10,'markerfacecolor','w')
hold off
% axis tight
ylim([min(nodes(:,2))-0.2,max(nodes(:,2))+0.2])
xlim([min(nodes(:,1))-0.2,max(nodes(:,1))+0.2])
axis square
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Graph')
simple_figure()
% ------------------------------------------------------------------------------

