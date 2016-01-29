function X=addframe(X)
X(1,:,:)=0;
X(end,:,:)=0;
X(:,1,:)=0;
X(:,end,:)=0;
end