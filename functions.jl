########
## 1D ##
########

heaviside(x)=0.5*(sign(x)+1)
pulse(t,ti,tf)=heaviside(t-ti)-heaviside(t-tf)

function pulse_train(t,tis,width)
    out=0.0
    for i=1:length(tis)
        out+=pulse(t,tis[i],tis[i]+width)
    end
    return out
end

S(x)=tanh(x)

function basic_1D!(du,u,p,t)
    
    τ = p[1]
    k = p[2]
    I = p[3]
    
    du[1]=1/τ*(-u[1]+S(k*u[1]+I(t)))
end

function σ_basic_1D!(du,u,p,t)
    
    du[1]=p[4]
end



###############
## 1D exc 1o ##
###############

function basic_1D_exc_1o!(du,u,p,t)
    
    τ = p[1]
    I = p[2]
    
    x = u[1]
    xs = u[2]

    k = p[4]
    ε = p[5]
    
    du[1]=1/τ*(-x+S(k*x+I(t)-xs))
    du[2] = ε*(x-xs)
end

function σ_basic_1D_exc_1o!(du,u,p,t)
    du[1]=p[3]
    du[2]=0.0
end


###########################
## 1D exc 1o, burst ##
###########################

function basic_1D_exc_1o_multi!(du,u,p,t)
    
    τ = p[1]
    I = p[2]
    
    x = u[1]
    xs = u[2]
    xu = u[3]

    k = p[4]
    ε = p[5]
    
    xs0 = p[6]
    ku = p[7]
    εu = p[8]
    
    du[1]=1/τ*(-x+S(k*x+I(t)-(xs-xs0)^2 -ku*xu))
    du[2] = ε*(x-xs)
    du[3] = εu*(x-xu)
end


function σ_basic_1D_exc_1o_multi!(du,u,p,t)
    du[1]=p[3]
    du[2]=0.0
    du[3]=0.0
end


###########################
## 1D exc 1o, burst, net ##
###########################

function basic_1D_exc_1o_multi_net!(du,u,p,t)
    
    τ = p[1]
    I = p[2]
    
    x = @view u[1:N]
    xs = @view u[N+1:2N]
    xu = @view u[2N+1:3N]
    xsyn = @view u[3N+1:4N]
    
    dx = @view du[1:N]
    dxs = @view du[N+1:2N]
    dxu = @view du[2N+1:3N]
    dxsyn = @view du[3N+1:4N]

    k = p[4]
    ε = p[5]
    
    xs0 = p[6]
    ku = p[7]
    εu = p[8]
    
    τsyn = p[9]
    
    dx[:]=1/τ*(-x+S.(k.*x + A*xsyn + I(t) -(xs.-xs0).^2 - ku.*xu))
    dxs[:] = ε*(x-xs)
    dxu[:] = εu*(x-xu)
    dxsyn[:] = 1/τsyn*(max.(0,x)-xsyn)
    
end


function σ_basic_1D_exc_1o_multi_net!(du,u,p,t)
    
    x = @view u[1:N]
    xs = @view u[N+1:2N]
    xu = @view u[2N+1:3N]
    xsyn = @view u[3N+1:4N]
    
    dx = @view du[1:N]
    dxs = @view du[N+1:2N]
    dxu = @view du[2N+1:3N]
    dxsyn = @view du[3N+1:4N]
    
    dx[:] .= p[3]
    dxs[:] .= 0.0
    dxu[:] .= 0.0
    dxsyn[:] .= 0.0
end



###############
## 1D exc 2o ##
###############

Sk(x,K,θ) = K*x^4/(θ^4+x^4)

function basic_1D_exc_2o!(du,u,p,t)
    
    τ = p[1]
    I = p[2]
    
    x = u[1]
    ks = u[2]

    k0 = p[4]
    Kx = p[5]
    K = p[6]
    ε = p[7]
    
    du[1]=1/τ*(-x+S((k0+Kx*x^2-ks)*x+I(t)))
    du[2] = ε*((K*x)^4-ks)
end

function σ_basic_1D_exc_2o!(du,u,p,t)
    du[1]=p[3]
    du[2]=0.0
end


#####################
## 1D exc 2o STIFF ##
#####################

function basic_1D_exc_2o_stiff!(du,u,p,t)
    
    τ = p[1]
    I = p[2]
    
    x = u[1]
    xs = u[2]

    k0 = p[4]
    Kx = p[5]
    ε = p[6]
    
    du[1]=1/τ*(-x+S((k0+Kx*x^2)*x+I(t)-xs))
    du[2] = ε*(x-xs)
end



###############
## 3D exc 3o ##
###############

function basic_3D_exc_3o!(du,u,p,t)
    
    x1 = u[1]
    x2 = u[2]
    x3 = u[3]
    ks1 = u[4]
    ks2 = u[5]
    ks3 = u[6]
    
    τ = p[1]
    I1 = p[2]
    I2 = p[3]
    I3 = p[4]

    k0 = p[6]
    Kx = p[7]
    K = p[8]
    ε = p[9]
    
    du[1]=1/τ*(-x1+S((k0+Kx*x1^2-ks1)*(-x2-x3)+I1(t)))
    du[2]=1/τ*(-x2+S((k0+Kx*x2^2-ks2)*(-x1-x3)+I2(t)))
    du[3]=1/τ*(-x3+S((k0+Kx*x3^2-ks3)*(-x1-x2)+I3(t)))
    
    du[4] = ε*((K*x1)^4-ks1)
    du[5] = ε*((K*x2)^4-ks2)
    du[6] = ε*((K*x3)^4-ks3)
end

function σ_basic_3D_exc_3o!(du,u,p,t)
    du[1]=p[5]
    du[2]=p[5]
    du[3]=p[5]
    
    du[4] = 0.0
    du[5] = 0.0
    du[6] = 0.0
end


function basic_3D_exc_3o_1Dfb!(du,u,p,t)
    
    x1 = u[1]
    x2 = u[2]
    x3 = u[3]
    ks = u[4]
    
    τ = p[1]
    I1 = p[2]
    I2 = p[3]
    I3 = p[4]

    k0 = p[6]
    Kx = p[7]
    K = p[8]
    ε = p[9]
    
    du[1]=1/τ*(-x1+S((k0+Kx*(x1^2+x2^2+x3^2)-ks)*(-x2-x3)+I1(t)))
    du[2]=1/τ*(-x2+S((k0+Kx*(x1^2+x2^2+x3^2)-ks)*(-x1-x3)+I2(t)))
    du[3]=1/τ*(-x3+S((k0+Kx*(x1^2+x2^2+x3^2)-ks)*(-x1-x2)+I3(t)))
    
    du[4] = ε*( K^4*(x1^4+x2^4+x3^4) - ks)
end


function σ_basic_3D_exc_3o_1Dfb!(du,u,p,t)
    du[1]=p[5]
    du[2]=p[5]
    du[3]=p[5]
    
    du[4] = 0.0
end





#####################
## Net 2D ##
#####################

f(v,k1,k2,I)=[-0*v[1]+tanh(k1*v[1]+I[1]);
                  -0*v[2]+tanh(k2*v[2]+I[2])]

g(x,k1,k2,I)=A*f(B*x,k1,k2,CC*I)

function one_slowD_dep!(du,u,p,t)
    
    k1 = p[1]
    k2 = p[2]
    
    I1 = p[3]
    I2 = p[4]
    
    I(t) = [I1(t);
            I2(t)]
    
    τ = p[7]
    
    x = [u[1];
         u[2]]
    
    du[1]=(-u[1]+g(x,k1,k2,I(t))[1]+0*I1(t))/τ
    du[2]=(-u[2]+g(x,k1,k2,I(t))[2]+0*I2(t))/τ
end

function σ_one_slowD_dep!(du,u,p,t)
    
    σ1 = p[5]
    σ2 = p[6]
    
    du[1]=σ1
    du[2]=σ2
end


#####################
## Net generic ##
#####################

fN(v,k,I)=tanh.(k.*v + I) #k,v,I are N-dimensional vectors

gN(x,k,I)=A*fN(B*x,k,CC*I) # A, B, C are N x N

function generic_net_2o!(du,u,p,t)
    
    k = p[1]
    
    I = p[2]
    
    τ = p[3]
    
    x = u
    
    du[:] = ( -x + gN(x,k,I(t))  )/τ
end

function σ_generic_net_2o!(du,u,p,t)
    
    σ = p[4]
    
    du[:] = σ*ones(N)
end


##################################
## Neural Field Kernel Building ##
##################################

function build_kernel_NF(σ,m,anys,periodic)
    A = zeros(m,m)
    for i=1:m
        for j=1:m
            if periodic
                A[i,j]=exp(-abs(sin(pi*mod(i-j,m)/m))/σ)
            else
                A[i,j]=exp(-abs(i-j)/(m*σ))*(1-anys*(j>i))
            end
        end
    end
    for i=1:m
        #A[i,i]=0
    end
    A=A/maximum(A)
    return A
end