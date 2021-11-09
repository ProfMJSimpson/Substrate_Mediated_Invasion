using Plots
pyplot()



gr()
default()
default(
    framestyle      =:axis,
    grid            =:y,
    axis            =:x,
    tick_direction  = :out,
    foreground_color_border = "#aaa",
    foreground_color_axis = "#aaa",
)
default(
    fontfamily      ="Helvetica",
    guidefontsize   = 12,
    titlefontsize   = 15,
)
colors = [["#fdbb84","#ef6548","#990000"],   # Blue, WM983b
          ["#99d8c9","#41ae76","#00441b"],   # Green, WM793b
          ["#cccccc","#969696","#525252"]]   # Grey, 451Lu

c_sph  = ["#A7D670","#F5473E","#B077C7"]


L=300.0
h=3.0
dt=0.01
N=Int(L/h)+1
xloc=zeros(N)
for i in 1:N
    xloc[i]=0.0+(i-1)*h
end
DD=400.0
K=1.0
r=0.6
r1=1.0/1.
r2=1.0/1.
tmax=14
maxsteps=Int(tmax/dt)
step1=Int(4/dt)
step2=Int(7/dt)
step3=Int(10/dt)
step4=Int(14/dt)
tt=0.0
u=K*ones(N,N)
pu=zeros(N,N)
us=zeros(N,N)
uplot=zeros(N,N,5)
splot=zeros(N,N,5)

s=r1*K*ones(N,N)/r2
smax=r1*K/r2


for i=2:N-1
    for j=2:N-1
        u[i,j]=0.0
        s[i,j]=0.0
    end
end




for k in 1:maxsteps
#println(k)


pu=u
ps=s

for i in 2:N-1
    for j in 2:N-1
u[i,j]=pu[i,j]+dt*(DD*((s[i,j]+s[i+1,j])*(u[i+1,j]-u[i,j])-(s[i,j]+s[i-1,j])*(u[i,j]-u[i-1,j])+(s[i,j]+s[i,j+1])*(u[i,j+1]-u[i,j])-(s[i,j]+s[i,j-1])*(u[i,j]-u[i,j-1]))/(2*h^2)+r*u[i,j]*(1-u[i,j]/K))
    end
end


for i in 1:N
    for j in 1:N
s[i,j]=ps[i,j]+dt*(r1*u[i,j]-r2*s[i,j])
    end
end

if k==1 
p1=contour(xloc,xloc,u,c=:winter,fill=true,aspect_ratio=1,xticks = [0,150,300],yticks = [0,150,300],xlabel="x",ylabel="y")
p2=contour(xloc,xloc,s,c=:spring,fill=true,aspect_ratio=1,xticks = [0,150,300],yticks = [0,150,300],xlabel="x",ylabel="y")
p3=plot(xloc,u[50,:],label=false,xticks = [0,150,300],yticks = [0,0.5,1.0],xlabel="x",ylabel="u(x,150,0)")
p4=plot(xloc,s[50,:],label=false,xticks = [0,150,300],yticks = [0,0.5,1.0],xlabel="x",ylabel="s(x,150,0)")
p5=plot(p1,p2,p3,p4,layout=(2,2))
display(p5)
savefig(p5,"Hole_Closing_t0.svg") 
elseif k==step1 
p1=contour(xloc,xloc,u,c=:winter,fill=true,aspect_ratio=1,xticks = [0,150,300],yticks = [0,150,300],xlabel="x",ylabel="y")
p2=contour(xloc,xloc,s,c=:spring,fill=true,aspect_ratio=1,xticks = [0,150,300],yticks = [0,150,300],xlabel="x",ylabel="y")
p3=plot(xloc,u[50,:],label=false,xticks = [0,150,300],yticks = [0,0.5,1.0],xlabel="x",ylabel="u(x,150,4)")
p4=plot(xloc,s[50,:],label=false,xticks = [0,150,300],yticks = [0,0.5,1.0],xlabel="x",ylabel="s(x,150,4)")
p5=plot(p1,p2,p3,p4,layout=(2,2))
display(p5)
savefig(p5,"Hole_Closing_t4.svg") 
elseif k==step2
p1=contour(xloc,xloc,u,c=:winter,fill=true,aspect_ratio=1,xticks = [0,150,300],yticks = [0,150,300],xlabel="x",ylabel="y")
p2=contour(xloc,xloc,s,c=:spring,fill=true,aspect_ratio=1,xticks = [0,150,300],yticks = [0,150,300],xlabel="x",ylabel="y")
p3=plot(xloc,u[50,:],label=false,xticks = [0,150,300],yticks = [0,0.5,1.0],xlabel="x",ylabel="u(x,150,7)")
p4=plot(xloc,s[50,:],label=false,xticks = [0,150,300],yticks = [0,0.5,1.0],xlabel="x",ylabel="s(x,150,7)")
p5=plot(p1,p2,p3,p4,layout=(2,2))
display(p5)
savefig(p5,"Hole_Closing_t7.svg")
elseif k==step3
p1=contour(xloc,xloc,u,c=:winter,fill=true,aspect_ratio=1,xticks = [0,150,300],yticks = [0,150,300],xlabel="x",ylabel="y")
p2=contour(xloc,xloc,s,c=:spring,fill=true,aspect_ratio=1,xticks = [0,150,300],yticks = [0,150,300],xlabel="x",ylabel="y")
p3=plot(xloc,u[50,:],label=false,xticks = [0,150,300],yticks = [0,0.5,1.0],xlabel="x",ylabel="u(x,150,10)")
p4=plot(xloc,s[50,:],label=false,xticks = [0,150,300],yticks = [0,0.5,1.0],xlabel="x",ylabel="s(x,150,10)")
p5=plot(p1,p2,p3,p4,layout=(2,2))
display(p5)
savefig(p5,"Hole_Closing_t10.svg")
elseif k==step4
p1=contour(xloc,xloc,u,c=:winter,fill=true,aspect_ratio=1,xticks = [0,150,300],yticks = [0,150,300],xlabel="x",ylabel="y")
p2=contour(xloc,xloc,s,c=:spring,fill=true,aspect_ratio=1,xticks = [0,150,300],yticks = [0,150,300],xlabel="x",ylabel="y")
p3=plot(xloc,u[50,:],label=false,xticks = [0,150,300],yticks = [0,0.5,1.0],xlabel="x",ylabel="u(x,150,14)")
p4=plot(xloc,s[50,:],label=false,xticks = [0,150,300],yticks = [0,0.5,1.0],xlabel="x",ylabel="s(x,150,14)")
p5=plot(p1,p2,p3,p4,layout=(2,2))
display(p5)
savefig(p5,"Hole_Closing_t14.svg")
end

end
