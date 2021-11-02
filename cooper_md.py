from utility import *
import matplotlib.pyplot as plt
from numpy import linspace
from scipy.integrate import odeint

class Cooper_md():
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
            "time": [0.0],
        }

        # propensity functions
        self.propensities = {
            # patient gets colonized
            0: lambda d: self.init_paras['r_hcw_p'] * d["p_s"][-1] * d["h_cln"][-1]
                         / self.init_paras['N_h'],

            # HCW gets colonized
            1: lambda d: self.init_paras['r_p_hcw']* d["h_s"][-1] * d["p_cln"][-1]
                         / self.init_paras['N_h'],

            # colonized HCW washes hands to eliminate pathogen
            2: lambda d: self.init_paras['r_hw'] * d["h_cln"][-1],

            # colonized patient is detected
            3: lambda d: self.init_paras['r_det'] * d["p_cln"][-1],

            # patient is regularly removed
            4: lambda d: self.init_paras['r_rm'] * d["p_cln"][-1],

            # Proportion of admissions already colonized
            5: lambda d: self.init_paras['pro_cln'] *
                (self.init_paras['r_rm'] *
                (d["p_s"][-1] + d["p_cln"][-1]+ self.init_paras['r_det'] * d["p_cln"][-1]) +
                self.init_paras['r_det'] * d["p_cln"][-1])
        }

        # change in species for each propensity
        self.stoichiometry = {
            0: {"p_s": -1, "p_cln": 1, "h_s": 0, "h_cln": 0},
            1: {"p_s": 0, "p_cln": 0, "h_s": -1, "h_cln": 1},
            2: {"p_s": 0, "p_cln": 0, "h_s": 1, "h_cln": -1},
            3: {"p_s": 1, "p_cln": -1, "h_s": 0, "h_cln": 0},
            4: {"p_s": 1, "p_cln": -1, "h_s": 0, "h_cln": 0},
            5: {"p_s": -1, "p_cln": 1, "h_s": 0, "h_cln": 0},
        }

    def __call__(self, plot = True,save = True):
        plt.figure(figsize=(10,10), dpi=500)
        fig = plt.gcf()
        fig.suptitle(self.title, fontsize=16)

        axes_pc = plt.subplot(211)
        axes_pc.set_ylabel("number of infected patients")
        axes_pc.set_xlabel('$t$')
        axes_pc.set_xlim([0, self.max_t])
        axes_pc.set_ylim([0, self.init_paras['N_p']])

        axes_hc = plt.subplot(212)
        axes_hc.set_ylabel("number of colonized HCWs")
        axes_hc.set_xlabel('$t$')
        axes_hc.set_xlim([0, self.max_t])
        axes_hc.set_ylim([0, self.init_paras['N_h']])

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
            axes_pc.step(trajectory["time"], trajectory["p_cln"], color="orange")
            axes_hc.step(trajectory["time"], trajectory["h_cln"], color="orange")
            trajectories += 1
            if trajectories == self.max_ite:
                break

        # numerical solution using an ordinary differential equation solversir
        t = linspace(0, self.max_t, num=self.max_t*10)
        y0 = (self.init_paras['p_s'], self.init_paras['p_cln'],
              self.init_paras['h_s'], self.init_paras['h_cln'])

        solution = odeint(self.differential_SIR, y0, t,
                          args=(self.init_paras['r_hcw_p'], self.init_paras['r_p_hcw'],
                                self.init_paras['r_hw'], self.init_paras['r_det'],
                                self.init_paras['r_rm'], self.init_paras['pro_cln']))

        solution = [[row[i] for row in solution] for i in range(4)]

        # plot numerical solution
        axes_pc.plot(t, solution[1], color="black")
        axes_hc.plot(t, solution[3], color="black")
        if plot:
            plt.show()
        if save:
            plt.savefig(self.title+'.png', dpi=800)

    def differential_SIR(self, n, t, r_hcw_p, r_p_hcw, r_hw, r_det, r_rm, pro_cln):
        dps_dt = -pro_cln * (r_rm * (n[0] + n[1]) + r_det * n[1]) \
                 - r_hcw_p * n[0] * n[3] / self.init_paras['N_p'] + (r_det+r_rm) * n[1]
        dpc_dt = pro_cln * (r_rm * (n[0] + n[1]) + r_det * n[1]) \
                 + r_hcw_p * n[0] * n[3] / self.init_paras['N_p'] - (r_det+r_rm) * n[1]
        dhs_dt = -r_p_hcw * n[2] * n[1] / self.init_paras['N_p'] + r_hw * n[3]
        dhc_dt = r_p_hcw * n[2] * n[1] / self.init_paras['N_p'] - r_hw * n[3]
        return dps_dt, dpc_dt, dhs_dt, dhc_dt