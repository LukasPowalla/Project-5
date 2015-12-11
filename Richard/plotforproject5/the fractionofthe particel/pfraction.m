fid=fopen('kinPotverlet50.txt')
    s1 = textscan(fid,'%f %f')
fclose(fid)

x1=s1{1};
y1=s1{2};

fid=fopen('kinPotverlet100.txt')
    s2 = textscan(fid,'%f %f')
fclose(fid)

x2=s2{1};
y2=s2{2};
fid=fopen('kinPotverlet200.txt')
    s3 = textscan(fid,'%f %f')
fclose(fid)

x3=s3{1};
y3=s3{2};
fid=fopen('kinPotverlet400.txt')
    s4 = textscan(fid,'%f %f')
fclose(fid)

x4=s4{1};
y4=s4{2};

fid=fopen('kinPotverlet800.txt')
    s5 = textscan(fid,'%f %f')
fclose(fid)

x5=s5{1};
y5=s5{2};
y2 = y2/1000;
y3 = y3/1000;
y4 = y4/1000;
hold on
%plot(y1,x1,'red')
plot(y2,x2,'blue')
plot(y3,x3,'green')
plot(y4,x4,'yellow')
%plot(y5,x5,'red')
legend('N=100','N=200','N=400');



xlabel('t/t_{crunch}')
ylabel('f^p(t)','Rotation',0)
%title('The fraction of the particles with positive energy as a funciton of time')