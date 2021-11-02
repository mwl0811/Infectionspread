from utility import *
import matplotlib.pyplot as plt
from numpy import linspace
from scipy.integrate import odeint

class Model():
    def __init__(
            self, init_paras, max_t, max_ite, title
    ):
        self.init_paras = init_paras
        self.max_t = max_t
        self.max_ite = max_ite
        self.title = title
        # initial species counts and sojourn times
        self.initital_conditions = {
            "p_s": [self.init_paras['p_s']],
            "p_cln": [self.init_paras['p_cln']],
            "h_s": [self.init_paras['h_s']],
            "h_cln": [self.init_paras['h_cln']],
            "nur_s": [self.init_paras['nur_s']],
            "nur_cln": [self.init_paras['nur_cln']],
            "time": [0.0],
        }


        # propensity functions
        self.propensities = {
            # patient gets colonized through contacts from HCWs
            0: lambda d: self.init_paras['r_hcw_p'] * d["p_s"][-1] * d["h_cln"][-1] / self.init_paras['N_h'],
            # patient gets colonized through contacts from nurses
            1: lambda d: self.init_paras['r_nur_p'] * d["p_s"][-1] * d["nur_cln"][-1] / self.init_paras['N_nur'],

            # HCW gets colonized
            2: lambda d: self.init_paras['r_p_hcw']* d["h_s"][-1] * d["p_cln"][-1] / self.init_paras['N_h'],

            # colonized HCW washes hands to eliminate pathogen
            3: lambda d: self.init_paras['r_hw_hcw'] * d["h_cln"][-1],

            # Nurse gets colonized
            4: lambda d: self.init_paras['r_p_nur']* d["nur_s"][-1] * d["p_cln"][-1] / self.init_paras['N_nur'],

            # colonized Nurse washes hands to eliminate pathogen
            5: lambda d: self.init_paras['r_hw_nur'] * d["nur_cln"][-1],

            # colonized patient is detected
            6: lambda d: self.init_paras['r_det'] * d["p_cln"][-1],

            # patient is regularly removed from the ward
            7: lambda d: self.init_paras['r_rm'] * d["p_cln"][-1],

            # Proportion of admissions already colonized
            8: lambda d: self.init_paras['pro_cln'] * (self.init_paras['r_rm'] * (d["p_s"][-1] + d["p_cln"][-1] +
                        self.init_paras['r_det'] * d["p_cln"][-1]) + self.init_paras['r_det'] * d["p_cln"][-1]),
        }

        # change in species for each propensity
        self.stoichiometry = {
            0: {"p_s": -1, "p_cln": 1,  "h_s": 0,  "h_cln": 0,  "nur_s":0,  "nur_cln":0 },
            1: {"p_s": -1, "p_cln": 1,  "h_s": 0,  "h_cln": 0,  "nur_s":0,  "nur_cln":0 },
            2: {"p_s": 0,  "p_cln": 0,  "h_s": -1, "h_cln": 1,  "nur_s":0,  "nur_cln":0 },
            3: {"p_s": 0,  "p_cln": 0,  "h_s": 1,  "h_cln": -1, "nur_s":0,  "nur_cln":0 },
            4: {"p_s": 0,  "p_cln": 0,  "h_s": 0,  "h_cln": 0,  "nur_s":-1, "nur_cln":1 },
            5: {"p_s": 0,  "p_cln": 0,  "h_s": 0,  "h_cln": 0,  "nur_s":1,  "nur_cln":-1},
            6: {"p_s": 1,  "p_cln": -1, "h_s": 0,  "h_cln": 0,  "nur_s":0,  "nur_cln":0 },
            7: {"p_s": 1,  "p_cln": -1, "h_s": 0,  "h_cln": 0,  "nur_s":0,  "nur_cln":0 },
            8: {"p_s": -1, "p_cln": 1,  "h_s": 0,  "h_cln": 0,  "nur_s":0,  "nur_cln":0 }
        }


    def __call__(self, plot = True,save = True):
        plt.figure(figsize=(10,10), dpi=500)
        fig = plt.gcf()
        fig.suptitle(self.title, fontsize=16)

        axes_p_cln = plt.subplot(311)
        axes_p_cln.set_ylabel("number of infected patients")
        axes_p_cln.set_xlabel('$t$')
        axes_p_cln.set_xlim([0, self.max_t])
        axes_p_cln.set_ylim([0, self.init_paras['N_p']])

        axes_hcw_cln = plt.subplot(312)
        axes_hcw_cln.set_ylabel("number of colonized HCWs")
        axes_hcw_cln.set_xlabel('$t$')
        axes_hcw_cln.set_xlim([0, self.max_t])
        axes_hcw_cln.set_ylim([0, self.init_paras['N_h']])

        axes_nur_cln = plt.subplot(313)
        axes_nur_cln.set_ylabel("number of colonized nurses")
        axes_nur_cln.set_xlabel('$t$')
        axes_nur_cln.set_xlim([0, self.max_t])
        axes_nur_cln.set_ylim([0, self.init_paras['N_nur']])

        # instantiate the epidemic SSA model container
        epidemic = SSAModel(
            self.initital_conditions,
            self.propensities,
            self.stoichiometry
        )

        # instantiate the SSA container with model
        epidemic_generator = SSA(epidemic, self.max_t)

        # simulate and plot 30 trajectories
        trajectories = 0
        for trajectory in epidemic_generator.direct():
            axes_p_cln.step(trajectory["time"], trajectory["p_cln"], color="orange")
            axes_hcw_cln.step(trajectory["time"], trajectory["h_cln"], color="orange")
            axes_nur_cln.step(trajectory["time"], trajectory["nur_cln"], color="orange")
            trajectories += 1
            if trajectories == self.max_ite:
                break

        # numerical solution using an ordinary differential equation solversir
        t = linspace(0, self.max_t, num=self.max_t*10)
        y0 = (self.init_paras['p_s'], self.init_paras['p_cln'],
              self.init_paras['h_s'], self.init_paras['h_cln'],
              self.init_paras['nur_s'], self.init_paras['nur_cln'])

        solution = odeint(self.differential_SIR, y0, t,
                          args=(self.init_paras['r_hcw_p'], self.init_paras['r_p_hcw'],
                                self.init_paras['r_hw_hcw'],
                                self.init_paras['r_nur_p'], self.init_paras['r_p_nur'],
                                self.init_paras['r_hw_nur'],
                                self.init_paras['r_det'],
                                self.init_paras['r_rm'], self.init_paras['pro_cln']))

        solution = [[row[i] for row in solution] for i in range(6)]

        # plot numerical solution
        axes_p_cln.plot(t, solution[1], color="black")
        axes_hcw_cln.plot(t, solution[3], color="black")
        axes_nur_cln.plot(t, solution[5], color="black")

        fig = plt.gcf()
        if plot:
            plt.show()
        if save:
            fig.savefig(self.title+'.png', dpi=800)

    def differential_SIR(self, n, t,  r_hcw_p, r_p_hcw, r_hw_hcw,  r_nur_p, r_p_nur, r_hw_nur, r_det, r_rm, pro_cln):
        dps_dt = -pro_cln * (r_rm * (n[0] + n[1]) + r_det * n[1]) \
                 - r_nur_p * n[0] * n[5] / self.init_paras['N_nur'] \
                 - r_hcw_p * n[0] * n[3] / self.init_paras['N_h'] + (r_det+r_rm) * n[1]

        dpc_dt = pro_cln * (r_rm * (n[0] + n[1]) + r_det * n[1]) \
                 + r_nur_p * n[0] * n[5] / self.init_paras['N_nur'] \
                 + r_hcw_p * n[0] * n[3] / self.init_paras['N_h'] - (r_det+r_rm) * n[1]

        dhs_dt = -r_p_hcw * n[2] * n[1] / self.init_paras['N_h'] + r_hw_hcw * n[3]
        dhc_dt = r_p_hcw * n[2] * n[1] / self.init_paras['N_h'] - r_hw_hcw * n[3]

        dnurs_dt = -r_p_nur * n[4] * n[1] / self.init_paras['N_nur'] + r_hw_nur * n[5]
        dnurc_dt = r_p_nur * n[4] * n[1] / self.init_paras['N_nur'] - r_hw_nur * n[5]

        return dps_dt, dpc_dt, dhs_dt, dhc_dt, dnurs_dt, dnurc_dt