########
## 1D ##
########

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



###############
## 1D exc 3o ##
###############

function basic_1D_exc_3o!(du,u,p,t)
    
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

function σ_basic_1D_exc_3o!(du,u,p,t)
    du[1]=p[5]
    du[2]=p[5]
    du[3]=p[5]
    
    du[4] = 0.0
    du[5] = 0.0
    du[6] = 0.0
end




#####################
## Net 2D 1 slow D ##
#####################

f(v,k1,k2,I)=[-0*v[1]+tanh(k1*v[1]+I[1]);
                  -0*v[2]+tanh(k2*v[2]+I[2])]

const A = 
[0.5  0.3
0.1  0.4]
#Matrix(1.0*LinearAlgebra.I,2,2)

const B =
#[0.0 -1.0
#-1.0 0.0]
[0.4  0.5
0.1  0.5]

pure_dom = eigvecs(A*B)[:,2]
pure_ndom = eigvecs(A*B)[:,1]
U=[pure_dom pure_ndom]
pure_dom_t = (inv(U)[1,:])

const CC =
#Matrix(1.0*LinearAlgebra.I,2,2)
[0.28  0.45
0.98  0.18]

Cpure_dom_t=(A*CC)'*pure_dom_t

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