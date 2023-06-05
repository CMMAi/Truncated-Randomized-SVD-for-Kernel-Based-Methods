function do_Fig05()
% This program compares the CPU time for rSVD and SVD 
%truncation are excluded
%%    
clear; 
close all
warning off
global RBFscale
global RBFpar
global RBFtype
RBFtype='g';
RBFpar=0.5;
RBFscale=0.1; %assume a scale parameter to construct the kernel matrix
%% inputs
n = round(logspace(log10(2000), log10(10000), 30));
n=n';
nt = 400; % number of test point
%% get the domain boundary
%initialization
toc_tsvd=zeros(length(n),1);
toc_trsvd=zeros(length(n),1);
%%%%%%%%%%%%%%%%%%%%%
%% LOOP
for j=1:length(n)
    disp(sprintf('n: %d',n(j)))
    % get the points
    [Pint, Pcntr, ~, ~, ~]=getPoints01(n(j),n(j),nt,0,1);
    A = kermat(Pint,Pcntr);
    
    tic
    [U,S,V] = svd(A);
    toc_tsvd(j) = toc;
    
    tic
    [Ur,Sr,Vr] = rsvd0(A,round(n(j)/2));
    toc_trsvd(j) = toc;
    
end
%%
%% plot
figure()
plot(n,toc_tsvd,'k-*','LineWidth',2);hold on
plot(n,toc_trsvd,'r-o','LineWidth',2)
xlabel('n')
ylabel('CPU time (s)')
legend('tSVD','trSVD','Location', 'northwest')
set(gca,'FontSize',16)
% title('CPU time (s)')3
str1=sprintf('Ex3_CPU_plot_tsvd_trsvd');
% saveas(gcf,str1,'eps') 
saveas(gcf,str1,'fig') 
