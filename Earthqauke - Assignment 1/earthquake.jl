using Random

function earthquake()
    # explosion detector (using spiral coorindate system) # run this file from julia by typing:]
    # julia> using PyPlot
    # julia> include("earthquake.jl")
    # julia> earthquake()
    # define the coordinate system:
    S=5000 # number of points on the spiral 
    x=zeros(S); y=zeros(S)
    for s=1:S
        theta=50*2*pi*s/S; r=s/S 
        x[s]=r*cos(theta); y[s]=r*sin(theta) 
    end
    plot(x,y,".")
    # define the locations of the detection stations on the surface
    # Also define what value on each sensor would be generated by an explostion at internal location s 
    N=10 # number of stations
    x_sensor=zeros(N); y_sensor=zeros(N)
    v=zeros(S,N)
    for sensor=1:N
        theta_sensor=2*pi*sensor/N
        x_sensor[sensor]=cos(theta_sensor); y_sensor[sensor]=sin(theta_sensor)
        for s=1:S
            v[s,sensor]=value(x[s],y[s],x_sensor[sensor],y_sensor[sensor]) # explosion value
        end
    end
    sd=1 # standard deviation of the Gaussian noise
    # Make the explosion data:
    true_s=rand(1:S) # true location of the explosion 
    theta=50*2*pi*true_s/S
    r=true_s/S
    x_true=r*cos(theta); y_true=r*sin(theta)
    # Get the noisy sensor values that will be observed for this explosion: 
    val=zeros(N)
    val_clean=zeros(N) # unknown clean values (just for interest)
    for sensor=1:N 
        val_clean[sensor]=value(x_true,y_true,x_sensor[sensor],y_sensor[sensor]) 
        val[sensor]=val_clean[sensor]+sd*randn()
    end
    figure()
    plot(1:N,val,label="noisy observed sensor measurements") 
    plot(1:N,val_clean,label="clean (unknown) sensor measurements") 
    legend()
    # Perform inference p(location|observed sensor values) given these sensor values
    logp=zeros(S)
    for s=1:S
        for sensor=1:N
            logp[s] += -0.5*(val[sensor]-v[s,sensor])^2/(sd^2) # Gaussian distribution 
        end
    end
    p=exp.(logp.-maximum(logp)) # do exponentiation (and avoid over/underflow) 
    p=p/sum(p) # normalise
    # plot the posterior and most likely location of the explosion:
    maxp,maxind =findmax(p)
    figure()
    for s=1:S
        plot(x[s],y[s],".",color=(1-(p[s]/maxp))*[1,1,1]) 
    end
    for theta=0:0.01:2*pi
        plot(cos(theta),sin(theta),".",color=[0,0,0])
    end
    for sensor=1:N
        plot(x_sensor[sensor],y_sensor[sensor],"o",color=[1,0,0])
    end
    plot(x_true,y_true,"rx",markersize=20,label="true") 
    plot(x[maxind],y[maxind],"m+",markersize=20,label="estimated (most likely)") 
    plot(sum(p.*x),sum(p.*y),"go",markersize=10,label="estimated (average)") 
    legend()
end
    
function value(x_true,y_true,x_sensor,y_sensor)
    return 1/(0.1+ (x_true-x_sensor)^2 + (y_true-y_sensor)^2) 
end