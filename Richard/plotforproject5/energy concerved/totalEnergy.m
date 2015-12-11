%fid=fopen('kinPotVerletE.txt')
fid=fopen('kinPotVerlet1.txt')
    s = textscan(fid,'%f %f %f')
fclose(fid)


x=s{1};
y=s{2};
x = x./x(1);
y = y/1000;

hold on
plot(y,x,'red')

xlabel('t/t_{crunch}')
ylabel('E/E_0')
%title('total energy over time with the smoothing function')