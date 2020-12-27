# parameters to control the order of the simulation
#param n_ens; #number of ens availables
#param e1; #number of experiments
#param e2; #required inputs in each experiment
param grid {(g1,g2) in {1 .. 6, 1 .. 9}}; # grid that controls the entire process

#param s_number; # number of scenarios
#param reps; # number of repetitions

# lenght of the time period since 2003
#param y0; #1: 2003, 8: 2010, 15: 2017, 22: 2024
#param y1; #1: 2003, 8: 2010, 15: 2017, 22: 2024

# -------- Sets --------
set bmi := {1 .. 20};
set age := {1 .. 7};
#set year := {1 .. 22} ordered;
set year := {22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1} ordered;
set ens := {1 .. 3} ordered;

# -------- Variables --------
var a {(i,j) in {bmi,age}} >= 0, <= 1; # Increase bmi
var b {(i,j) in {bmi,age}} >= 0, <= 1; # Decrease bmi
var c {(i,j) in {bmi,age}} >= 0, <= 1; # Remain bmi
var U {(i,j,k) in {bmi,age,year}} integer, >= 0; # Estimated population
var BETA_AUM {(i,j,k) in {bmi,age,year}} >= 0, <= 1; # Increase and change age range
var BETA_DIS {(i,j,k) in {bmi,age,year}} >= 0, <= 1; # Decrease and change age range
var BETA_MAN {(i,j,k) in {bmi,age,year}} >= 0, <= 1; # Remain and change age range
var GAMA_AUM {(i,j,k) in {bmi,age,year}} >= 0, <= 1; # Increase and remain age range
var GAMA_DIS {(i,j,k) in {bmi,age,year}} >= 0, <= 1; # Decrease and remain age range
var GAMA_MAN {(i,j,k) in {bmi,age,year}} >= 0, <= 1; # Remain and remain age range

# -------- Parameters --------
# Demographic
param pop {(i,j,d) in {bmi,age,ens}}; # known populations from ENS
param y {(j,k) in {age,year}}; # Proportion of population at age range boundaries
param m {(i,j) in {bmi,age}}; # Population growth
param w {j in age}; # Age range weights

# Feasible bounds
param amin {(i,j) in {bmi,age}}; # Increase bmi, lower bound
param amax {(i,j) in {bmi,age}}; # Increase bmi, upper bound
param bmin {(i,j) in {bmi,age}}; # Decrease bmi, lower bound
param bmax {(i,j) in {bmi,age}}; # Decrease bmi, upper bound
param cmin {(i,j) in {bmi,age}}; # Remain bmi, lower bound
param cmax {(i,j) in {bmi,age}}; # Remain bmi, upper bound

# Auxiliar, to measure the time
param t0;

# Parameters to control the population to work with
param p0 {(i,j) in {bmi,age}}; # initial population (known)
param p1 {(i,j) in {bmi,age}}; # objective population (known)
#param p2 {(i,j) in {bmi,age}}; # objective population (estimated)

param ord_year_init; # initial year
param ord_ens_init;
param ord_year_obj; # objective year
param ord_ens_obj;

# Funcion objetivo que considera periodo 2003-2010 (p1-p2)
#minimize z : sum {(i,j) in {bmi,age}} ((T[i,j]-U[i,j,8])^2)*w[j];
minimize z : sum {(i,j,k) in {bmi,age,year}: ord(k) == ord_year_obj} ((p1[i,j]-U[i,j,k])^2)*w[j];

#------------------------------Restricciones--------------------------

/* # initial population (known)
let {(i,j,d) in {bmi,age,ens}: ord(d) == ord_ens_init} p0[i,j] := pop[i,j,d];
# objective population (known)
let {(i,j,d) in {bmi,age,ens}: ord(d) == ord_ens_obj} p1[i,j] := pop[i,j,d]; */

# initial population (known)
subject to res_init_pop {(i,j,k) in {bmi,age,year}: ord(k) == ord_year_init}:
U[i,j,k] = p0[i,j];

/* # objective population (estimated)
subject to res_f_obj_pop {(i,j,k) in {bmi,age,year}: ord(k) == ord_year_obj}:
U[i,j,k] = p2[i,j]; */

#-------Esquina superior izquierda
subject to restriccion1 {(i,j,k) in {bmi,age,year}: i==1 and j==1 and k!=last(year)}:
((1-y[j,k])*(U[i,j,k]*(1-a[i,j])+b[i,j]*U[i+1,j,k]))*(1-m[i,j])=U[i,j,next(k,year)];

#--------Esquina inferior izquierda
subject to restriccion3 {(i,j,k) in {bmi,age,year}: i==card(bmi) and j==1 and k!=last(year)}:
((1-y[j,k])*(U[i,j,k]*(1-b[i,j])+a[i,j]*U[i-1,j,k]))*(1-m[i-1,j])=U[i,j,next(k,year)];

#--------Esquina superior derecha
subject to restriccion5 {(i,j,k) in {bmi,age,year}: i==1 and j==card(age) and k!=last(year)}:
(U[i,j,k]*(1-a[i,j])+b[i,j]*U[i+1,j,k]+(1-b[i,j-1]-a[i,j-1])*y[j-1,k]*U[i,j-1,k]+b[i,j-1]*y[j-1,k]*U[i+1,j-1,k])*(1-m[i,j])=U[i,j,next(k,year)];

#---------Esquina inferior derecha
subject to restriccion7 {(i,j,k) in {bmi,age,year}: i==card(bmi) and j==card(age) and k!=last(year)}:
(U[i,j,k]*(1-b[i,j])+a[i,j]*U[i-1,j,k]+(1-b[i,j-1]-a[i,j-1])*y[j-1,k]*U[i,j-1,k]+a[i,j-1]*y[j-1,k]*U[i-1,j-1,k])*(1-m[i,j])=U[i,j,next(k,year)];

#----------Borde superior
subject to restriccion9 {(i,j,k) in {bmi,age,year}: i==1 and j>1 and j<card(age) and k!=last(year)}:
((1-y[j,k])*(U[i,j,k]*(1-a[i,j])+b[i,j]*U[i+1,j,k])+y[j-1,k]*((1-a[i,j-1])*U[i,j-1,k]+b[i,j-1]*U[i+1,j-1,k]))*(1-m[i,j])=U[i,j,next(k,year)];

#----------Borde izquierdo
subject to restriccion11 {(i,j,k) in {bmi,age,year}: i>1 and i<card(bmi) and j==1 and k!=last(year)}:
((1-y[j,k])*(U[i,j,k]*(1-b[i,j]-a[i,j])+a[i,j]*U[i-1,j,k]+b[i,j]*U[i+1,j,k]))*(1-m[i,j])=U[i,j,next(k,year)];

#----------Borde inferior
subject to restriccion13 {(i,j,k) in {bmi,age,year}: i==card(bmi) and j>1 and j<card(age) and k!=last(year)}:
((1-y[j,k])*(U[i,j,k]*(1-b[i,j])+a[i,j]*U[i-1,j,k])+y[j-1,k]*((1-b[i,j-1])*U[i,j-1,k]+a[i,j-1]*U[i-1,j-1,k]))*(1-m[i,j])=U[i,j,next(k,year)];

#----------Borde derecho
subject to restriccion15 {(i,j,k) in {bmi,age,year}: i>1 and i<card(bmi) and j==card(age) and k!=last(year)}:
(U[i,j,k]*(1-b[i,j]-a[i,j])+a[i,j]*U[i-1,j,k]+b[i,j]*U[i+1,j,k]+y[j-1,k]*(a[i,j-1]*U[i-1,j-1,k]+(1-b[i,j-1]-a[i,j-1])*U[i,j-1,k]+b[i,j-1]*U[i+1,j-1,k]))*(1-m[i,j])=U[i,j,next(k,year)];

#----------Centro
subject to restriccion17 {(i,j,k) in {bmi,age,year}: i>1 and i<card(bmi) and j>1 and j<card(age) and k!=last(year)}:
((1-y[j,k])*(U[i,j,k]*(1-b[i,j]-a[i,j])+a[i,j]*U[i-1,j,k]+b[i,j]*U[i+1,j,k])+y[j-1,k]*(a[i,j-1]*U[i-1,j-1,k]+(1-b[i,j-1]-a[i,j-1])*U[i,j-1,k]+b[i,j-1]*U[i+1,j-1,k]))*(1-m[i,j])=U[i,j,next(k,year)];

#------------------------------Restricciones de tasas------------------------------
subject to restriccion21 {i in bmi}:
sum {j in age} (a[i,j] + b[i,j] + c[i,j]) = 1;

subject to restriccion22 {(i,j) in {bmi,age}}:
a[i,j] >= amin[i,j];

subject to restriccion23 {(i,j) in {bmi,age}}:
a[i,j] <= amax[i,j];

subject to restriccion24 {(i,j) in {bmi,age}}:
b[i,j] >= bmin[i,j];

subject to restriccion25 {(i,j) in {bmi,age}}:
b[i,j] <= bmax[i,j];

subject to restriccion26 {(i,j) in {bmi,age}}:
c[i,j] >= cmin[i,j];

subject to restriccion27 {(i,j) in {bmi,age}}:
c[i,j] <= cmax[i,j];

subject to restriccion28 {(i,j,k) in {bmi,age,year}}:
BETA_DIS[i,j,k] = y[j,k]*b[i,j];

subject to restriccion29 {(i,j,k) in {bmi,age,year}}:
BETA_MAN[i,j,k] = y[j,k]*c[i,j];

subject to restriccion30 {(i,j,k) in {bmi,age,year}}:
BETA_AUM[i,j,k] = y[j,k]*a[i,j];

subject to restriccion31 {(i,j,k) in {bmi,age,year}}:
GAMA_DIS[i,j,k] = (1-y[j,k])*b[i,j];

subject to restriccion32 {(i,j,k) in {bmi,age,year}}:
GAMA_MAN[i,j,k] = (1-y[j,k])*c[i,j];

subject to restriccion33 {(i,j,k) in {bmi,age,year}}:
GAMA_AUM[i,j,k] = (1-y[j,k])*a[i,j];
