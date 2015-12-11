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



hold on
y2 = y2 /1000;
%plot(y1,x1,'red')
plot(y2,x2,'blue')
plot(y2,x3,'green')
plot(y2,x4,'yellow')
%plot(y5,x5,'red')


legend('N=200','N=100','N=400');



xlabel('t/t_{crunch}')
ylabel('f^p(t)','Rotation',0)
%title('The fraction of the particles with positive energy as a funciton of time without smoothing function')