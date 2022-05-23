using Pkg
Pkg.add("HomotopyContinuation")
Pkg.add("DelimitedFiles")
Pkg.add("LinearAlgebra")
Pkg.add("SymPy")


using HomotopyContinuation
using DelimitedFiles
using LinearAlgebra


@var p11 p12 p13  p14 p15 p16 p21 p22 p23 p24 p25 p26 p31 p32 p33 p34 p35 p36 

P1=[p11 p12 p13; p12 p14 p15; p13 p15 p16]
P2=[p21 p22 p23; p22 p24 p25; p23 p25 p26]
P3=[p31 p32 p33; p32 p34 p35; p33 p35 p36]


### Example 1 parameters 
#
A=[1 0 0;0 1 0;0 0 1]; B1=[1;0;0]; B2=[0;1;0]; B3=[0;0;1]
R1=1; R2=1; R3=1
Q1=[1 -1 0;-1 1 0;0 0 0]; Q2=[1 -1 0;-1 2 -1;0 -1 1]; Q3=[0 0 0;0 1 -1;0 -1 1] 
S1=B1*R1*B1'; S2=B2*R2*B2'; S3=B3*R3*B3'
A1=A-S2*P2-S3*P3; A2=A-S1*P1-S3*P3; A3=A-S1*P1-S2*P2
stblcntindx=Int16(1)
#
are1=A1'*P1+P1*A1+Q1-P1*S1*P1
are2=A2'*P2+P2*A2+Q2-P2*S2*P2
are3=A3'*P3+P3*A3+Q3-P3*S3*P3
#

### Coupled Algebriac Riccati Equations  (18 polynomial equations in 18 variables)
###
equations = [are1[1,1]; are1[1,2]; are1[1,3]; are1[2,2]; are1[2,3]; are1[3,3];
+are2[1,1]; are2[1,2]; are2[1,3]; are2[2,2]; are2[2,3]; are2[3,3];
+are3[1,1]; are3[1,2]; are3[1,3]; are3[2,2]; are3[2,3]; are3[3,3]]

###  Solve the polynomails equaitons using Homotopy Continuation method 
###
F = System(equations)
allsolutions = solve(F)

### Real solutions of CARE
###
sol=real_solutions(allsolutions)
writedlm( "allsolutions.csv",  sol, ',')
sz=size(sol)
numberofsolutions=sz[1]

###  To find the number of stabilizing solutions of CARE
###
using SymPy

p11,p12,p13,p14,p15,p16 = symbols("p11 p12 p13 p14 p15 p16")
p21,p22,p23,p24,p25,p26 = symbols("p21 p22 p23 p24 p25 p26")
p31,p32,p33,p34,p35,p36 = symbols("p31 p32 p33 p34 p35 p36")
#
P1=[p11 p12 p13; p12 p14 p15; p13 p15 p16]
P2=[p21 p22 p23; p22 p24 p25; p23 p25 p26]
P3=[p31 p32 p33; p32 p34 p35; p33 p35 p36]
#
stblindx=zeros(Int16,numberofsolutions)
for k in 1:numberofsolutions
    sol1=sol[k]
    sP1=float(P1.subs([(p11,sol1[1]),(p12,sol1[2]),(p13,sol1[3]),(p14,sol1[4]),(p15,sol1[5]),(p16,sol1[6])]))
    sP2=float(P2.subs([(p21,sol1[7]),(p22,sol1[8]),(p23,sol1[9]),(p24,sol1[10]),(p25,sol1[11]),(p26,sol1[12])]))
    sP3=float(P3.subs([(p31,sol1[13]),(p32,sol1[14]),(p33,sol1[15]),(p34,sol1[16]),(p35,sol1[17]),(p36,sol1[18])]))
    Acl=float(A-S1*sP1-S2*sP2-S3*sP3)
    ev=real(eigvals(Acl))
    flag=Float64(all(<(0),ev))    
    if flag == 1
        stblindx[stblcntindx]=k
        global stblcntindx=stblcntindx+1
    end
end
numberofStablesolutions = stblcntindx-1
writedlm( "stablesolutions.csv",  sol[stblindx[1:stblcntindx-1]], ',')

sol1=sol[stblindx[1]]
sP1=float(P1.subs([(p11,sol1[1]),(p12,sol1[2]),(p13,sol1[3]),(p14,sol1[4]),(p15,sol1[5]),(p16,sol1[6])]))
sP2=float(P2.subs([(p21,sol1[7]),(p22,sol1[8]),(p23,sol1[9]),(p24,sol1[10]),(p25,sol1[11]),(p26,sol1[12])]))
sP3=float(P3.subs([(p31,sol1[13]),(p32,sol1[14]),(p33,sol1[15]),(p34,sol1[16]),(p35,sol1[17]),(p36,sol1[18])]))

### To find the number of solutions of CARE which sastisfy the structural constriants (7)
### 
### The CARE solutions satisfy the structural constriants if p13=0 and p33=0
concntindx =1
conindx=zeros(Int16,numberofsolutions)
for k in 1:numberofsolutions
    sol1=sol[k]
    sp13=sol1[3]
    sp33=sol1[15]
    flag = (sp13 == 0) && (sp33 == 0) 
    if flag == 1
        conindx[concntindx]=k
        global concntindx=concntindx+1
    end
end

numberofsolutionsConstraints = concntindx-1


println("\n\n Number of solutions of CARE -- ", numberofsolutions, "\n")
println("\n solutions are written in the file allsolutions.csv\n")
println("\n\n Number of stablizing solutions of CARE -- ", numberofStablesolutions, "\n")
println("\n Stabilizing solution of CARE is ")
println("\n P1 -- ", sP1, "\n P2 -- ", sP2, "\n P3 --", sP3)
println("\n stable solutions are written in the file stablesolutions.csv\n")
println("\n\n Number of stablizing solutions of CARE satisfying the structural constraints (7) -- ",numberofsolutionsConstraints, "\n")
