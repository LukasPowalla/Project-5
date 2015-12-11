fid=fopen('kinPotVerletE.txt')
    s = textscan(fid,'%f %f %f')
fclose(fid)

x=s{1};
y=s{2};
z=s{3};
z=z/z(1);

x=x*(-1);
x = x +100;

hold on
plot(x,z,'.');


xlabel('particles ejected form the system')
ylabel('E/E_0','Rotation',90)
%title('enery of the system vs particles ejected form the system with N = 100')