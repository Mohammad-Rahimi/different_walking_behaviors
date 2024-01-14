%% Cs are based on reconstruction using non negative least square matlab:
source_files=["08000normalsecond",
"08000wide2second",
"08000widersecond",
 "08015normal4second",
"08015widesecond",
"08015widersecond",
 "08030normal4second",
"08030wide4second",
"08030wider3second",
 "11000normalsecond",
"11000widesecond",
"11000widersecond",
 "11015normal3second",
"11015widesecond",
"11015widersecond",
 "11030normal7second",
"11030wide2second",
"11030wider5second",
"14500normalsecond",
"14500widesecond",
"14500widersecond",
"14515normalsecond",
"14515widesecond",
"14515wider2second",
"14530normalsecond",
"14530widesecond",
"14530wider2second"
];
time=[0
0.01
0.02
0.03
0.04
0.05
0.06
0.07
0.08
0.09
0.1
0.11
0.12
0.13
0.14
0.15
0.16
0.17
0.18
0.19
0.2
0.21
0.22
0.23
0.24
0.25
0.26
0.27
0.28
0.29
0.3
0.31
0.32
0.33
0.34
0.35
0.36
0.37
0.38
0.39
0.4
0.41
0.42
0.43
0.44
0.45
0.46
0.47
0.48
0.49
0.5
0.51
0.52
0.53
0.54
0.55
0.56
0.57
0.58
0.59
0.6
0.61
0.62
0.63
0.64
0.65
0.66
0.67
0.68
0.69
0.7
0.71
0.72
0.73
0.74
0.75
0.76
0.77
0.78
0.79
0.8
0.81
0.82
0.83
0.84
0.85
0.86
0.87
0.88
0.89
0.9
0.91
0.92
0.93
0.94
0.95
0.96
0.97
0.98
0.99
1];
theta=time*2*pi;
data2=importdata("compareCs_reconstfinalJan.xls")
data22=data2.data
excel(1,1)="nnlsqr"
syndata1=importdata('11000normalsecond_19-Jul-2023.mat');
act1=importdata(strcat("11000normalsecond",".sto"))
act2=importdata(strcat(filename,".sto"))
cr=syndata1.synergyOutRight(4).C;
cl=syndata1.synergyOutLeft(4).C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:2
    
    if i==1
      data1=importdata("all_recons_C_r_finalJan.xls");
      side="Right"
      refC=cr
      a=2
    else if i==2
      data1=importdata("all_recons_C_l_finalJan.xls");
      side="Left";
      refC=cl
      a=1
        end
    end
    for n=1:length(source_files)
        for j=1:4
            xx=refC(j,:)
            yy=data1((n-1)*4+j,:)
            avg(j)=mean(data1((n-1)*4+j,:));
            maxactiv(j)=max(data1((n-1)*4+j,:));
            excel(n+1,1)=source_files(n);
            excel(n+1,j+1)=avg(j);
            excel(1,j+1)="average";
            excel(n+1,j+5)=maxactiv(j);
            excel(1,j+5)="max";
            excel(1,j+9)="R"
            excel(n+1,j+9)=corr(xx',yy');
            
            excel(1,j+13)="CoA"
            yyCoa=(atan2((yy*sin(theta)),(yy*cos(theta))));
            if yyCoa<0
               yyCoa=yyCoa+2*pi;
            end
            yyCoa=yyCoa/(2*pi);
            excel(n+1,j+13)=yyCoa;
            excel(1,j+17)="dwt"
            excel(n+1,j+17)=dtw(xx',yy');
            excel(1,j+21)="MI"
            excel(n+1,j+21)=mi_cont_cont(xx',yy');
        end
    end
    xlswrite(strcat("CsMeanMaxRCoA_finaljan",side),excel)  
    
end
%-------------------------------------------------------------------------------
% mi_cont_cont: mutual information using the first method in [1]
%
% Syntax: mi = mi_cont_cont(x, y, k)
%
% Inputs: 
%     x - input vector (length-N)
%     y - input vector (length-N)
%     k - number of nearest neighbours (default 3)
%
% Outputs: 
%     mi - mutual information estimate (scalar)
%
% Example:
%     % set the number of nearest neighbours:
%     k = 3;
%     
%     % generate low-frequency random signals:
%     x = randn(1, 1000);
%     y = randn(1, 1000);
%     
%     % calculate MI between continuous x and discrete y:
%     mi = mi_cont_cont(x, y, k);
%     fprintf('MI = %g\n', mi);
%     
% 
% Requires:
%     'knnsearch' from the statistics toolbox
% 
%
%  [1] Kraskov, A., Stögbauer, H., & Grassberger, P. (2004). Estimating mutual
%  information. Physical Review E, 69(6), 16. https://doi.org/10.1103/PhysRevE.69.066138
% John M. O' Toole, University College Cork
% Started: 05-08-2020
%
% last update: Time-stamp: <2020-08-12 13:34:38 (otoolej)>
%-------------------------------------------------------------------------------
function mi = mi_cont_cont(x, y, k)
if(nargin < 3 || isempty(k)), k = 5; end
if(length(x) ~= length(y))
    error('x and y should be the same length');
end
N = length(x);
x = x(:);
y = y(:);
xy = [x y];
%---------------------------------------------------------------------
% 1. calculate the k-th nearest neighbour:
%---------------------------------------------------------------------
nn_mdlx = createns(xy, 'Distance','chebychev');
[~, dist_kth] = knnsearch(nn_mdlx, xy, 'k', k + 1);
dist_kth = dist_kth(:, k + 1);
%---------------------------------------------------------------------
%  2. find all points within distance:
%  find all the points (all classes) within the k-th distance :
%---------------------------------------------------------------------
% sort vectors first:
[xs, idx] = sort(x);
[ys, idy] = sort(y);
[~, idx] = sort(idx);
[~, idy] = sort(idy);
nx_range = rangesearch_var_radius(xs, idx, dist_kth);
ny_range = rangesearch_var_radius(ys, idy, dist_kth);
%---------------------------------------------------------------------
% 3. calcuate MI (see [1])
%---------------------------------------------------------------------
mi = psi(k) + psi(N) - mean(psi(nx_range + 1) + psi(ny_range + 1));
% correct if <0
mi(mi < 0) = 0;
end
%-------------------------------------------------------------------------------
% rangesearch_var_radius: 'rangesearch' for varying distances (radii)
%
% Syntax: nx_range = rangesearch_var_radius(x, idx, d)
%
% Inputs: 
%     x    - vector of data points (length-N array)
%     idx  - sequence of indices of the data points from x (length-N array)
%     d    - vector of distances (length-N array)
%
% Outputs: 
%     nx_range - number of samples within distance d(n) for x(idx(n)), 
%                for n=1,2,...,N (length-N array)
%
% Example:
%     % define dummy variables:
%     N = 1000;
%     x = randn(1, N);
%     d = abs(randn(1, N)) + 0.1;
%     
%     % return the number of samples within the radius d:
%     nx_range = rangesearch_var_radius(x, 1:N, d);
%     
%
% John M. O' Toole, University College Cork
% Started: 09-08-2020
%
% last update: Time-stamp: <2020-08-12 15:33:42 (otoolej)>
%-------------------------------------------------------------------------------
function nx_range = rangesearch_var_radius(x, idx, d)
N = length(x);
    
nx_range = zeros(1, N);
for n = 1:N
    nx_range(n) = rangesearch_1point(x, idx(n), d(n), N);
end
end   
function p = rangesearch_1point(x, ix0, d, N)
%---------------------------------------------------------------------
% how many points within (less than) the distance d for point x(ix0)
%---------------------------------------------------------------------
x0 = x(ix0);
p = 0;
for n = (ix0 + 1):N
    if(abs(x(n) - x0) >= d)
        break;
    else
        p = p + 1;
    end
end
for n = (ix0 - 1):-1:1
    if(abs(x(n) - x0) >= d)
        break;
    else
        p = p + 1;
    end
end
end

