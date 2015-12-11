fid=fopen('kinPotVerlet.txt')
    s = textscan(fid,'%f %f %f %f')
fclose(fid)

x=s{1};
y=s{2};
z=s{3};
r=s{4};
x=x/1000;
y = y./z(1);
z = z./z(1);
r=r./z(1);

hold on
plot(x,y,'yellow')
plot(x,z,'blue')
%plot(x,r,'green')

xlabel('t/t_{crunch}')
ylabel('E/E_0','Rotation',90)
title('The behvier og kinetic and potential energy with 100 objects')
legend('kinetic energy','potential')