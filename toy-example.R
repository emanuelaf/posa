# Toy example
source('codice/functions.R')
N = 10
pi_init = rep(0.2, N)
y = c(0,1,0,1,0,0,0,1,1,0)
cond = 0
pi_i = calculate_pi_i(pi_init = pi_init, y=y, cond=cond)
#pi_ik = calculate_pi_ik(pi_init = pi_init, y = y, cond = cond, pi_i = pi_i)
posa_ht_est = numeric(100000)
posa_pht_est = numeric(100000)
var_pht_posa = calculate_var_pht_posa(pi_init=pi_init, pi_i=pi_i, y=y, cond=cond, length(y))[[1]]
var_ht_posa = calculate_var_ht_posa(pi_i = pi_i, pi_init = pi_init, 
                                    y = y, N = length(y), cond = cond)

for (i in 1:100000) {
   s = posa(pik = pi_init, Y=y, cond=cond)
   posa_ht_est[i] = t_posa_ht(s.pi_i = pi_i[s==1], s.y = y[s==1])/length(y)
   posa_pht_est[i] = t_posa_pht(s = s, y.s = y[s==1], pik = pi_init, cond=cond)/length(y)
}
(mean(posa_ht_est)-sum(y)/length(y))/(sum(y)/length(y))*100
(mean(posa_pht_est)-sum(y)/length(y))/(sum(y)/length(y))*100
(var_pht_posa-var(posa_pht_est))/(var(posa_pht_est))*100
(var_ht_posa-var(posa_ht_est))/(var(posa_ht_est))*100

