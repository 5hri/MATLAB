function m = cpv(eigvalue,getpercent)
% cumulative percentage variance 
a=sum(eigvalue);
[Row,Col]=size(eigvalue);
n=max(Row,Col);

eigvalue = sort(eigvalue,'descend');
s=0;
for i=1:n
    s  =s + eigvalue(i)/a;
    if s > getpercent
        m = i;
        break;
    end
end   