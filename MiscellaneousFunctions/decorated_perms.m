function G = decorated_perms(n,varargin)
% D=decorated_perms(N) returns a matrix, each row of which is a permutation
% of the numbers 1:N wherein fixed elements take one of two orientations
% for each fixed point (loop): clockwise (postive) or counterclockwise
% (negative).
% 
% D=decorated_perms(N,K) returns only those decorated permutations for
% which the sum of the number of exceedances and the number of
% counterclockwise loops is equal to K.  In this case, the number of rows
% of D is the number of totally positive cells in the positive Grassmannian
% Gr^+(n,k).  See Williams, L.K., 2005. Enumeration of totally positive
% Grassmann cells. Advances in Mathematics, 190(2), pp.319-342.
% 

P=flipud(perms(1:n));
nfact=length(P);
D=[];
for i=1:nfact
    afixed=find(P(i,:)==1:n);
    Dnew = P(i,:);
    if length(afixed)>0
        for k=1:length(afixed)
            Dnewpm=Dnew;
            Dnewpm(:,afixed(k))=-Dnewpm(:,afixed(k));
            Dnew = [ Dnew ; Dnewpm ];
        end
    end
    D = [ D ; Dnew ];
end

if nargin==1
    G=D;
else
    G=[];
    ktarget=varargin{1};
    for i=1:length(D)
        k=0;
        for j=1:n
            if abs(D(i,j))>j % exceedance  
                k=k+1;
            end
            if D(i,j)==-j % counter-clockwise loop
                k=k+1;
            end
        end
        if k==ktarget
            G = [ G ; D(i,:) ];
        end
    end
end

return





