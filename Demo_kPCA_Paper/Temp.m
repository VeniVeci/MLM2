[a,index]=sort(RE_P(:,end))

new=X_star(index',:)


[a,index]=sort(X_star(:,end));
new=RE_P(index',:)

plot(new(:,end))
% set(gca,'yscale','log');