fid=fopen('kinPotVerlet200.txt')
    s = textscan(fid,'%f %f %f %f')
fclose(fid)

x=s{1};
y=s{2};
z=s{3};
k=s{4};
z=z/z(1);
x=x/1000;

fid=fopen('kinPotVerlet100.txt')
    s2 = textscan(fid,'%f %f %f %f')
fclose(fid)

x2=s2{1};
y2=s2{2};
z2=s2{3};
k2=s2{4};
z2=z2/z2(1);
x2=x2/1000;


hold on
plot(x,z,'blue');
plot(x2,z2,'red');

xlabel('particles ejected form the system')
ylabel('E/E_0','Rotation',90)

xlabel('t/t_{crunch}')
legend('N=200','N=100')
%title('enery of the system vs particles ejected form the system with N = 100')