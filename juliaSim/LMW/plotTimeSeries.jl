include("TimeSeriesLMW.jl")
using Plots, Main.MySystem
using LinearAlgebra: norm
plotly()

upper_bound = 1e6
fp_threshold = 1e-2
periodic_threshold = 1e-2
time_threshold = 1e-2

## Decreasing upper_bound makes more of the plot cream
## Increasing fp_threshold makes more of the plot orange
## Increasing periodic_threshold makes more of the plot purple
## Increasing time_ threshold makes more of the plot purple

function checkFP(time_series)
    second_components = [row[2] for row in time_series]  # Extract the second component
    max_second_component = maximum(second_components)
    min_second_component = minimum(second_components)
    return max_second_component - min_second_component
end

function find_maxima(time_series, f,size)
    maxima_data = []
    derivatives = [row[size+2:end] for row in time_series]
    n = length(derivatives)

    for i in 1:n-1
        if derivatives[i][1] > 0 && derivatives[i+1][1] <= 0
            y = time_series[i][2:size+1]
            fy = time_series[i][size+2:end]
            t1 = time_series[i][1]
            t2 = time_series[i+1][1]
            d1 = derivatives[i][1]
            d2 = derivatives[i+1][1]
            te = (t2 - t1) * (0 - d1) / (d2 - d1)
            ye,fye,_ = rk(y,f,fy,te)

            if fye[1] > 0
                te = te + (te-t2) * (0-fye[1]) / (fye[1]-d2)
                ye, fye = rk(y,f,fy,te)
            else
                te = te*(0-d1)/(fye[1]-d1)
                ye, fye, _ = rk(y,f,fy,te)                
            end

            push!(maxima_data,[te,ye...])
        end
    end

    return maxima_data
end



function checkPeriodic(maxima_data, periodic_threshold, time_threshold, size)
    n = length(maxima_data)
    if n==0
        return 0
    end
    eps = periodic_threshold
    y1 = maxima_data[1][2:size+1]
    t1 = maxima_data[1][1]
    flag1 = 0
    flag2 = 0
    ind=0
    t2=0
    t3=0

    for i = 2:n
        if (norm(maxima_data[i][2:size+1] - y1) < eps)
            t2 = maxima_data[i][1]
            ind = i
            flag1 = 1
            break
        end
    end
    
    if flag1==0
        return 0
    end

    for i = ind+1:n
        if norm(maxima_data[i][2:size+1]-y1) < eps
            t3 = maxima_data[i][1]
            flag2=1
            break
        end
    end

    if flag2==0
        return 0
    end

    d1 = t2-t1
    d2 = t3-t2

    if abs(d1-d2) < time_threshold
        return 1
    end
    return 0
end

nd = 150
M = zeros(nd,nd)
X = range(-3, stop=3, length=nd)
Y = range(-3, stop=3, length=nd)
constant_params = (-1.0, 1.0, 1.0, -1.0, -1.8, -1.8)
#constant_params = (-1.0,1.0)


size=8



for (i,x) = enumerate(X)
    for (j,y) = enumerate(Y)
        params = (constant_params...,x,y)
        series, flag = integrate(rhs(params...),size; upper_bound=upper_bound)
        println(length(series))

        if flag ==1
            M[j,i] = 3
            continue
        end
        
        a = Int(floor(0.7*length(series)))
        diff = checkFP(series[a:end])
        
        if diff<fp_threshold
            M[j,i] = 2
            continue
        end
        maxi = find_maxima(series[a:end],rhs(params...),size)

        if checkPeriodic(maxi, periodic_threshold, time_threshold,size)==1
            M[j,i] = 1
            continue
        end
    end
end

#for (i,x) = enumerate(X)
#    for (j,y) = enumerate(Y)
#        if (x+y <= 0)&& (M[j,i] ==3)
#            M[j,i] =0
#        end
#    end
#end

parameterPlot = heatmap(X,Y,M, xlabel="j1",ylabel="j2", size=(600,600))
display(parameterPlot)
#abhiat07@gmail.com


# Parameter plotting