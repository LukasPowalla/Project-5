fid=fopen('radialDensity100.txt')
    s100 = textscan(fid,'%f %f')
fclose(fid)

x100=s100{1};
y100=s100{2};
x100 = x100./(100^(-1/3));
y100 = y100./100/100;


fid=fopen('radialDensity200.txt')
    s200 = textscan(fid,'%f %f')
fclose(fid)

x200=s200{1};
y200=s200{2};
x200 = x200./(200^(-1/3));
y200 = y200./200/200;



fid=fopen('radialDensity50.txt')
    s50 = textscan(fid,'%f %f')
fclose(fid)

x50=s50{1};
y50=s50{2};
x50 = x50./(50^(-1/3));
y50 = y50./50/50;



fid=fopen('radialDensity400.txt')
    s400 = textscan(fid,'%f %f')
fclose(fid)

x400=s400{1};
y400=s400{2};
x400 = x400./(400^(-1/3));
y400 = y400./400/400;


fid=fopen('radialDensity800.txt')
    s800 = textscan(fid,'%f %f')
fclose(fid)

x800=s800{1};
y800=s800{2};
x800 = x800./(800^(-1/3));
y800 = y800./800/800;



fid=fopen('radialDensity1000.txt')
    s1000 = textscan(fid,'%f %f')
fclose(fid)

x1000=s1000{1};
y1000=s1000{2};
x1000 = x1000./(1000^(-1/3));
y1000 = y1000./1000/1000;


x = logspace(-2,2);
y =  1./ (1+ (x./0.95).^4);


hold on ;
figure;
loglog(x,y,'black',x1000,y1000,'blue.',x800,y800,'red.',x400,y400,'yellow.',x100,y100,'green.')
%loglog(x2,y2,'blue.')
%loglog(x,y,'red.')
%plot(x200,y200,'black.')
%plot(x100,y100,'red.')
%plot(x50,y50,'blue.')


%loglog(x,y,'black');


%plot(x400,y400,'yellow.');
%loglog(x800,y800,'green.');
%loglog(x1000,y1000,'black.');


xlabel('r/N^{-1/3}')
ylabel('n(r)/N^2')

legend('1/(1+(r/0.95)^4)','N=1000','N=800','N=400','N=100');
