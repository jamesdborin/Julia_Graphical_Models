using Random
using PyPlot
using Printf


function value(x_true,y_true,x_sensor,y_sensor)
    return 1/(0.1+ (x_true-x_sensor)^2 + (y_true-y_sensor)^2) 
end


function readearthquakedata(filename)
    lines = open(filename) do file
        readlines(file)
    end
    [parse(Float64, i) for i in lines]
end


function spiralcoordinates(S, rate)
    x=zeros(S); y=zeros(S)
    for s=1:S
        theta=rate*2*pi*s/S;  r=s/S
        x[s]=r*cos(theta); y[s]=r*sin(theta)
    end
    return x, y
end


function sensorcoordinates(N)
    x_sensor=zeros(N); y_sensor=zeros(N)
    for sensor=1:N
        theta_sensor=2*pi*sensor/N
        x_sensor[sensor]=cos(theta_sensor); y_sensor[sensor]=sin(theta_sensor)
    end
    return x_sensor, y_sensor
end


function sensorexplosionvalues(S, N, x, y, x_sensor, y_sensor)
    v=zeros(S,N)
    for sensor=1:N
        for s=1:S
            v[s,sensor]=value(x[s],y[s],x_sensor[sensor],y_sensor[sensor]) # explosion value for some value function
        end
    end
    return v
end


function logprobabilities_twoexplosions(S, N, sd, val, v)
    logp = zeros(S,S)
    for s1=1:S
        for s2=1:S
            for sensor=1:N
                logp[s1,s2] += -0.5*(val[sensor]-v[s1,sensor]-v[s2,sensor])^2/(sd^2)
            end
        end
    end
    return logp
end


function logprobabilities_oneexplosion(S, N, sd, val, v)
    logp = zeros(S)
    for s=1:S
        for sensor=1:N
            logp[s] += -0.5*(val[sensor]-v[s,sensor])^2/(sd^2)
        end
    end
    return logp
end


function plotboundary()
    for theta=0:0.01:2*pi
        plot(cos(theta),sin(theta),".",color=[0,0,0])
    end
end


function plotsensors(x_sensor, y_sensor)
    for sensor=1:length(x_sensor)
        plot(x_sensor[sensor],y_sensor[sensor],"o",color=[1,0,0])
    end
end

function plotprobabilities(p, x, y)
    maxp, maxpind = findmax(p)
    for s=1:length(x)
        plot(x[s],y[s],".",color=(1-(p[s]/maxp))*[1,1,1]) 
    end
end


function plotposterior(filename, valuescale)
    S = 2000
    rate = 25
    sd = 0.2

    x, y = spiralcoordinates(S, rate)

    N=30 # number of stations

    x_sensor, y_sensor = sensorcoordinates(N)

    v = valuescale.*sensorexplosionvalues(S, N, x, y, x_sensor, y_sensor)

    val = readearthquakedata(filename)

    logp = logprobabilities_twoexplosions(S, N, sd, val, v)
    p=exp.(logp.-maximum(logp))
    p=p/sum(p)

    maxp, maxind = findmax(p)
    # Plot most likely locations
    figure()
    p1 = sum(p, dims=2)
    plotprobabilities(p1, x, y)
    plotboundary()
    plotsensors(x_sensor, y_sensor)
    plot(x[maxind[1]],y[maxind[1]],"m+",markersize=20,label="explosion 1") 
    plot(x[maxind[2]],y[maxind[2]],"m+",markersize=20,label="explosion 2") 

end

function plotposterioroneexplosion(filename, valuescale)
    S = 2000
    rate = 25
    sd = 0.2

    x, y = spiralcoordinates(S, rate)

    N=30 # number of stations

    x_sensor, y_sensor = sensorcoordinates(N)

    v = valuescale.*sensorexplosionvalues(S, N, x, y, x_sensor, y_sensor)

    val = readearthquakedata(filename)

    logp = logprobabilities_oneexplosion(S, N, sd, val, v)
    p=exp.(logp.-maximum(logp))
    p=p/sum(p)

    maxp, maxind = findmax(p)
    # Plot most likely locations
    figure()
    plotprobabilities(p, x, y)
    plotboundary()
    plotsensors(x_sensor, y_sensor)
    plot(x[maxind],y[maxind],"m+",markersize=20,label="explosion")
end

function hypothesislogratio(filename, valuescale)
    S = 2000
    rate = 25
    sd = 0.2

    x, y = spiralcoordinates(S, rate)

    N=30 # number of stations

    x_sensor, y_sensor = sensorcoordinates(N)

    v = valuescale.*sensorexplosionvalues(S, N, x, y, x_sensor, y_sensor)

    val = readearthquakedata(filename)

    logp2 = logprobabilities_twoexplosions(S, N, sd, val, v)
    maxp2 = maximum(logp2)
    p2=exp.(logp2.-maxp2)
    lp2 = log(sum(p2))+maxp2

    logp1 = logprobabilities_oneexplosion(S, N, sd, val, v)
    maxp1 = maximum(logp1)
    p1 = exp.(logp1.-maxp1)
    lp1 = log(sum(p1))+maxp1

    println(lp2)
    println(lp1)
    return lp2-lp1
end

function assignment1()
    firstdatafn = "EarthquakeExerciseData.txt"
    plotposterior(firstdatafn, 1.0)

    meandatafn = "EarthquakeExerciseMeanData.txt"
    plotposterior(meandatafn, 0.5)

    plotposterioroneexplosion(firstdatafn, 1.0)

    ratio = hypothesislogratio(firstdatafn, 1.0)
    println("Hypothesis log ratio for first data:")
    println("$(@sprintf("%.2f", ratio))")

    ratio2 = hypothesislogratio(meandatafn, 0.5)
    println("Hypothesis log ratio for mean data:")
    println("$(@sprintf("%.2f", ratio2))")
end
