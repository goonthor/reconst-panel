function data = scaleinput(d, fsize)
range = max(d)-min(d);
factor = max(abs(d)).*min(2*fsize./range, 1);
data = d./repmat(factor,numel(d(:,1)),1);