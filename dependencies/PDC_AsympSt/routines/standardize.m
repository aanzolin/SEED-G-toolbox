function c=standardize(x)

[m,n]=size(x);
if m > n, x=x.'; [m,n]=size(x); end;
for i=1:m
     c(:,i)=x(:,i)/std(x(:,i));
end