%Times plot

N2 = csvread('N2_alpha077_beta022_times.csv');
N6 = csvread('N6_alpha049_beta039_times.csv');
MCS = csvread('MCS.csv');
MCS = MCS(1,:);

%colormap
cc = hsv(6);

%N2plot
subplot(2,1,1);
hold all
for i = 1:1:6
    plot(log10(MCS),log10(N2(i,:))-log10(N2(1,:)),'o-')
end
xlabel('log10(MCS)');
ylabel('log10(Time[s]) - log10(BF, NLE)')
legend('(BF,NLE)','(BF,ALE)','(IS,NLE,NQF)','(IS,NLE,AQF)','(IS,ALE,NQF)','(IS,ALE,AQF)')

%N6plot
subplot(2,1,2);
hold all
for i = 1:1:6
    plot(log10(MCS),log10(N6(i,:))-log10(N6(1,:)),'o-')
end
xlabel('log10(MCS)');
ylabel('log10(Time[s]) - log10(BF, NLE)')