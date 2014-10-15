%monte carlo
counter = 1;
theta = linspace(0,2*pi,401);
c = cos(theta);
s = sin(theta);
hit = 0.0;
total = 0.0;
for k=1:1:1111
    subplot(1,3,counter)
    plot(c,s,'k-','LineWidth',2)
    hold on
    x = (rand - 0.5)*2;
    y = (rand - 0.5)*2;
    if (x^2 + y^2 < 1)
        hit = hit + 1;
    end
    total = total + 1;
    plot([x],[y],'x')
    if (k==10 | k==110) 
        counter = counter +1;
        subplot(1,3,counter)
        plot(c,s,'k-','LineWidth',2)
        pi_approx = 4.0*hit/total
        error = abs(pi_approx-pi)
        hold on
        hit = 0.0;
        total = 0.0;
    end
end
pi_approx = 4.0*hit/total
error = abs(pi_approx-pi)