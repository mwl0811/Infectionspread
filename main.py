from refined_md import Model
# number of HCWs
N_h = 10
# number of patients
N_p = 40
# number of nurses
N_nur = 10

# initial number of colonized patients
p_cln = 0
# initial number of susceptible patients
p_s = N_p - p_cln
# initial number of colonized HCWs
h_cln = 0
# initial number of susceptible HCWs
h_s = N_h - h_cln
# initial number of colonized nurses
nur_cln = 0
# initial number of susceptible nurses
nur_s = N_nur - nur_cln

# transmission rate =
# contact rate x probability of colonization at each contact
# assume probability of colonization at each contact = 0.1
# each patient requires 1 contacts from HCWs per day
# and requires 5 contacts from nurses per day


# HCW-patient transmission rate
r_hcw_p = 0.1
# patient-HCW transmission rate
r_p_hcw = 0.1
# hand washing rate of HCWs
r_hw_hcw = 6

# nur-patient transmission rate
r_nur_p = 0.5
# patient-nur transmission rate
r_p_nur = 0.5
# hand washing rate of nurses
r_hw_nur = 6

# detection rate of colonized patients
r_det = 0.1
# removal rate of colonized patients
r_rm = 0.1
# Proportion of admissions already colonized
pro_cln = 0.01



title = 'Simulation result'
# max time
max_t = 80
# max iteration
max_ite = 4

parameters = {
    'N_h': N_h, 'N_p': N_p, 'N_nur': N_nur,
    'p_cln': p_cln,'p_s': p_s,
    'h_cln': h_cln, 'h_s': h_s,
    'nur_cln': nur_cln, 'nur_s': nur_s,
    'r_hw_hcw': r_hw_hcw, 'r_hcw_p': r_hcw_p, 'r_p_hcw': r_p_hcw,
    'r_hw_nur': r_hw_nur, 'r_nur_p': r_nur_p, 'r_p_nur': r_p_nur,
    'r_det': r_det, 'r_rm': r_rm, 'pro_cln': pro_cln
}


model = Model(parameters, max_t, max_ite, title)
model(plot=True, save=True)
