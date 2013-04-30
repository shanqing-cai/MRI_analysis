function [PVAL,ADJ,NULL]=nbs_bct_sc(X, Y, testName, THRESH, K, TAIL)
%NBS_BST        Network-based statistic, as described in [1]. 
%
%     PVAL = NBS(X,Y,THRESH) performs the NBS for populations X and Y for a
%     T-statistic threshold of THRESH. The third dimension of X and Y 
%     references a particular member of the populations. The first two 
%     dimensions reference the connectivity value of a particular edge 
%     comprising the connectivity matrix. For example, X(i,j,k) stores the 
%     connectivity value corresponding to the edge between i and j for the
%     kth memeber of the population. PVAL is a vector of corrected p-values 
%     for each component identified. If at least one of the p-values is 
%     less than 0.05, then the omnibus null hypothesis can be rejected at 
%     5% significance. The null hypothesis is that the value of 
%     connectivity at each edge comes from distributions of equal mean 
%     between the two populations.
%
%     [PVAL,ADJ] = NBS(X,Y,THRESH) also returns an adjacency matrix 
%     identifying the edges comprising each component. Edges corresponding 
%     to the first p-value stored in the vector PVAL are assigned the value
%     1 in the adjacency matrix ADJ, edges corresponding to the second 
%     p-value are assigned the value 2, etc. 
%
%     [PVAL,ADJ,NULL] = NBS(X,Y,THRESH) also returns a vector of K samples 
%     from the the null distribution of maximal component size. 
%   
%     [PVAL,ADJ] = NBS(X,Y,THRESH,K) enables specification of the number of
%     permutations to be generated to estimate the empirical null
%     distribution of maximal component size. Default: K=1000. 
%
%     [PVAL,ADJ] = NBS(X,Y,THRESH,K,TAIL) enables specification of the type
%     of alternative hypothesis to test. If TAIL:
%     'both'  - alternative hypothesis is means are not equal (default)
%     'left'  - mean of population X < mean of population Y
%     'right' - mean of population X > mean of population Y 
%
%     ALGORITHM DESCRIPTION 
%     The NBS is a nonparametric statistical test used to isolate the 
%     components of an N x N undirected connectivity matrix that differ 
%     significantly between two distinct populations. Each element of the 
%     connectivity matrix stores a connectivity value and each member of 
%     the two populations possesses a distinct connectivity matrix. A 
%     component of a connectivity matrix is defined as a set of 
%     interconnected edges. 
%
%     The NBS is essentially a procedure to control the family-wise error 
%     rate, in the weak sense, when the null hypothesis is tested 
%     independently at each of the N(N-1)/2 edges comprising the 
%     connectivity matrix. The NBS can provide greater statistical power 
%     than conventional procedures for controlling the family-wise error 
%     rate, such as the false discovery rate, if the set of edges at which
%     the null hypothesis is rejected constitues a large component or
%     components.
%     The NBS comprises fours steps:
%     1. Perform a two-sample T-test at each edge indepedently to test the
%        hypothesis that the value of connectivity between the two
%        populations come from distributions with equal means. 
%     2. Threshold the T-statistic available at each edge to form a set of
%        suprathreshold edges. 
%     3. Identify any components in the adjacency matrix defined by the set
%        of suprathreshold edges. These are referred to as observed 
%        components. Compute the size of each observed component 
%        identified; that is, the number of edges it comprises. 
%     4. Repeat K times steps 1-3, each time randomly permuting members of
%        the two populations and storing the size of the largest component 
%        identified for each permuation. This yields an empirical estimate
%        of the null distribution of maximal component size. A corrected 
%        p-value for each observed component is then calculated using this
%        null distribution.
%
%     [1] Zalesky A, Fornito A, Bullmore ET (2010) Network-based statistic:
%         Identifying differences in brain networks. NeuroImage.
%         10.1016/j.neuroimage.2010.06.041
%
%     Written by: Andrew Zalesky, azalesky@unimelb.edu.au
%     Revised by: (2013-04-29) Shanqing Cai, shanqing.cai@gmail.com
%     




%Error checking
if nargin<3
    error('Not enough inputs\n');
end
if nargin<4
    K=1000;  
end
if nargin<5
    TAIL='both';
end    


%% -- Make sure that testName is valid -- %
if ~isequal(testName, 'ttest2') && ~isequal(testName, 'ranksum') ...
   && ~isequal(testName, 'lincorr')
    error('Unrecognized testName: %s', testName);
end

%%
[Ix,Jx,nx]=size(X);
if isequal(testName, 'ttest2') || isequal(testName, 'ranksum')    
    [Iy,Jy,ny]=size(Y);
    
    if any([Ix~=Jx,Iy~=Jy,Ix~=Jy])
        error('Matrices are not square, or are not of equal dimensions\n');
    end
else % -- lincorr -- %
    dimsY = length(size(Y));
    if dimsY ~= 2
        error('Wrong dimension in Y: %d', dimsY');        
    end
    
    Y = Y(:);
    [ny, wy] = size(Y);
    if wy ~= 1
        error('Under lincorr, Y must be a vector (not a matrix)');
    end
    if nx ~= ny;
        error('Length of Y (%d) does not match the number of subjects in X (%d)', ...
              ny, nx);
    end
end
     
    

%Number of nodes
N=Ix;

%Only consider elements above upper diagonal due to symmetry
ind=find(triu(ones(N,N),1));

%Number of edges
M=length(ind);

%Look up table
ind2ij=zeros(M,2);
[ind2ij(:,1),ind2ij(:,2)]=ind2sub([N,N],ind);

%Vectorize connectivity matrices
%Not necessary, but may speed up indexing
%Uses more memory since cmat temporarily replicates X
cmat=zeros(M,nx); 
for i=1:nx
    tmp=squeeze(X(:,:,i));
    cmat(:,i)=tmp(ind)';
end
clear X

if isequal(testName, 'ttest2') || isequal(testName, 'ranksum')
    pmat=zeros(M,ny); 
    for i=1:ny
        tmp=squeeze(Y(:,:,i));
        pmat(:,i)=tmp(ind)';
    end
else % -- lincorr -- %
    pmat = repmat(Y', length(ind), 1);
end
clear Y

%%

%Perform T-test at each edge
stat=zeros(M,1);
for i=1:M
    %[a,p_val,c,tmp]=ttest2(cmat(i,:),pmat(i,:));
    if isequal(testName, 'ttest2')
        tmp = ttest2_stat_only(cmat(i,:),pmat(i,:)); 
    elseif isequal(testName, 'ranksum')
        tmp = ranksum(cmat(i, :), pmat(i, :));
        tmp = sign(median(cmat(i, :)) - median(pmat(i, :))) * -log10(tmp);
    elseif isequal(testName, 'lincorr')
        [kk, r2, tmp] = lincorr(cmat(i, :), pmat(i, :));
        
        tmp = sign(kk(2)) * -log10(tmp);
    end
    stat(i)=tmp;
end

if strcmp(TAIL,'both')
    stat=abs(stat);
elseif strcmp(TAIL,'left')
    stat=-stat;
elseif strcmp(TAIL,'right')
else
    error('Tail option not recognized\n');
end

%Threshold 
ind_t=find(stat > THRESH);

%Suprathreshold adjacency matrix
ADJ=spalloc(N,N,length(ind_t)); 
ADJ(ind(ind_t))=1; 
ADJ=ADJ+ADJ'; 

bgl=0;
if exist('components')==2
    %Use components.m provided by MatlabBGL
    bgl=1;
end
%Find network components    
if bgl==1
    [a,sz]=components(ADJ);
else
    [a,sz]=get_components(ADJ);
end

%Convert size from number of nodes to number of edges
%Only consider components comprising more than one nodes (equivalent to at
%least one edge)
ind_sz=find(sz>1);
sz_links=zeros(1,length(ind_sz));
for i=1:length(ind_sz)
    nodes=find(ind_sz(i)==a); 
    sz_links(i)=sum(sum(ADJ(nodes,nodes)))/2;
    ADJ(nodes,nodes)=ADJ(nodes,nodes)*(i+1); 
end

%Subtract 1 to delete edges not comprising a component
%While 1 is also subtracted from edges comprising a component, this is 
%compensated by the (i+1) above.
ADJ(find(ADJ))=ADJ(find(ADJ))-1;

if ~isempty(sz_links)
    max_sz=max(sz_links);
else
    max_sz=0;
end
fprintf('Max component size is: %d\n',max_sz); 

%Empirically estimate null distribution of maximum compoent size by
%generating K independent permutations. 
fprintf('Estimating null distribution with permutation testing\n');
hit=0;
for k=1:K
    %Randomise
    if isequal(testName, 'ttest2') || isequal(testName, 'ranksum')
        d=zeros(M,nx+ny);
        d=[cmat,pmat];
        indperm=randperm(nx+ny); 
        d=d(:,indperm);
    else
        d = cmat;
        indperm = randperm(nx);
        d=d(:, indperm);
        d = [d, pmat];
    end
    
    %Perform T-test at each edge
    t_stat_perm=zeros(M,1);
    for i=1:M
        %[z1,z2,z3,tmp]=ttest2(d(i,1:nx),d(i,nx+1:nx+ny));
        if isequal(testName, 'ttest2')
            tmp = ttest2_stat_only(d(i,1:nx), d(i,nx+1:nx+ny));
        elseif isequal(testName, 'ranksum')
            tmp = ranksum(d(i,1:nx), d(i,nx+1:nx+ny));
            tmp = sign(median(d(i,1:nx)) - median(d(i,nx+1:nx+ny))) * -log10(tmp);
        elseif isequal(testName, 'lincorr')
            [kk, r2, tmp] = lincorr(d(i, 1 : nx), d(i, nx + 1 : nx + ny));
            tmp = sign(kk(2)) * -log10(tmp);
        end
            
        t_stat_perm(i)=tmp;
    end
    if strcmp(TAIL,'both')
        t_stat_perm=abs(t_stat_perm);
    elseif strcmp(TAIL,'left')
        t_stat_perm=-t_stat_perm;
    elseif strcmp(TAIL,'right')
    else
        error('Tail option not recognized\n');
    end
    
    %Threshold
    ind_t=find(t_stat_perm>THRESH);
    
    %Suprathreshold adjacency matrix
    adj_perm=spalloc(N,N,length(ind_t));
    adj_perm(ind(ind_t))=1;
    adj_perm=adj_perm+adj_perm';

    %Find size of network components   
    if bgl==1
        [a,sz]=components(adj_perm);
    else
        [a,sz]=get_components(adj_perm);
    end
    
    %Convert size from number of nodes to number of links 
    ind_sz=find(sz>1);
    sz_links_perm=zeros(1,length(ind_sz));
    for i=1:length(ind_sz)
        nodes=find(ind_sz(i)==a);
        sz_links_perm(i)=sum(sum(adj_perm(nodes,nodes)))/2;
    end
  
    if ~isempty(sz_links_perm)
        NULL(k)=max(sz_links_perm);
    else
        NULL(k)=0; 
    end
    if NULL(k)>=max_sz
        hit=hit+1;
    end
    fprintf('Perm %d of %d. Perm max is: %d. Observed max is: %d. P-val estimate is: %0.3f\n',k,K,NULL(k),max_sz,hit/k)
end

%Calculate p-values
for i=1:length(sz_links)
    PVAL(i)=length(find(NULL>=sz_links(i)))/K;
end




function t=ttest2_stat_only(x,y)

t=mean(x)-mean(y); 
n1=length(x); 
n2=length(y); 
s=sqrt(((n1-1)*var(x)+(n2-1)*var(y))/(n1+n2-2)); 
t=t/(s*sqrt(1/n1+1/n2)); 


