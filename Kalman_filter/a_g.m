theta = [0,5,10,15,20,25,30,35,40,45];

figure(1)
for i = 1:length(theta)
    subplot(2,5,i), loglog(logspace(-3,2,1000),1+(4*theta(i)*logspace(-3,2,1000))+(logspace(-3,2,1000).^2))
end

