speed = zeros(10, 1);
for i = 1:5;
    tic;
    A = rand(i*1000, i*1000, 'gpuArray');
    B = svd(A);
    B = gather(B);
    time1 = toc;
    
    tic; 
    C = rand(i*1000, i*1000);
    D = svd(C);
    time2 = toc;
    
    disp(time2/time1);
    speed(i) = time2/time1;
end
plot(speed);
    
    
    