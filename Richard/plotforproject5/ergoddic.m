fid=fopen('ergodic.txt')
    s = textscan(fid,'%f %f')
fclose(fid)

x=s{1};
y=s{2};

plot(y,x,'red')