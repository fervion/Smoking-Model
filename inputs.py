"""
Handles reading and manipulating inputs for the smoking model

@author: Taige Hou(thou1@partners.org)
"""

import xlrd, re
import numpy as np
from enums import *
import collections
######################################################################
#data structure for reading inputs

#tabs in the excel sheet. (sheetname, varname)
TABS = [("Simulation","sim"),
        ("Initialization", "init"),
        ("NatHist","nathist"),
        ("Smoking","smoking"),
        ("Pre-Events","pre_events"),
        ("Events","events"),
        ("Prophs","prophs"),
        ("Intervention", "intervention"),
        ("Costs","costs"),
        ("QoL","qol"),
        ]

#cells in excel sheet
#format is (txt label,varname,cell, celltype)
#inputname is name in the txt file
#if var is part of an array put varname in tuple (name, shape, index)
INPUTS = {}
INPUTS['sim'] = [("runsize", "runsize","D3", int),
                 ("disc_rate_year", "disc_rate_year","D4", float),
                 ("rand_seed", "rand_seed","E6",bool),
                 ("event_names", "event_names","B12;B{0}".format(NUM_EVENTS+11), str),
                 ("maxage","maxage","B19",int),
                 ]
def add_prev(sheet):
    prev = np.ndarray(((NUM_EVENTS, len(GENDERS), len(SS),len(SI),len(AGE_BRACKETS))))
    for event in range(NUM_EVENTS):
        for gender in GENDERS:
            offset = event*20+gender*10
            for si in SI:
                prev[event][gender][SS_NEVER][si] = xl_get_range(sheet, "Y{0}:AM{0}".format(7+offset), float)
            prev[event][gender][SS_FORMER] = xl_get_range(sheet, "Y{0}:AM{1}".format(8+offset, 10+offset),float)
            prev[event][gender][SS_CURRENT] = xl_get_range(sheet, "Y{0}:AM{1}".format(11+offset, 13+offset),float)

    return prev

def add_prev_load_text(textfile, inputs):
    prev = np.ndarray(((NUM_EVENTS, len(GENDERS), len(SS),len(SI),len(AGE_BRACKETS))))
    for event in range(NUM_EVENTS):
        for gender in GENDERS:
            offset = event*20+gender*10
            label = "prev_event\t{}\t{}".format(inputs.sim.event_names[event], SS_STRS[SS_NEVER])
            for si in SI:
                prev[event][gender][SS_NEVER][si] = text_get_range(textfile, "Y{0}:AM{0}".format(7+offset), label,float)
            label = "prev_event\t{}\t{}".format(inputs.sim.event_names[event], SS_STRS[SS_FORMER])
            prev[event][gender][SS_FORMER] = text_get_range(textfile, "Y{0}:AM{1}".format(8+offset, 10+offset),label,float)
            label = "prev_event\t{}\t{}".format(inputs.sim.event_names[event], SS_STRS[SS_CURRENT])
            prev[event][gender][SS_CURRENT] = text_get_range(textfile, "Y{0}:AM{1}".format(11+offset, 13+offset),label,float)

    return prev

def add_prev_write_text(textfile, inputs):
    prev=inputs.init.event_prev
    for event in range(NUM_EVENTS):
        for gender in GENDERS:
            for ss in SS:
                label = "prev_event\t{}\t{}\t".format(inputs.sim.event_names[event],SS_STRS[ss])
                textfile.write(label)
                if ss==SS_NEVER:
                    textfile.write("\t".join([str(_) for _ in prev[event,gender,ss,SI_LIGHT]]))
                else:
                    for si in SI:
                        if si!=0:
                            textfile.write("\n")
                        textfile.write("\t".join([str(_) for _ in prev[event,gender,ss,si]]))
                textfile.write("\n")


INPUTS['init'] = [("initage","age_dist","C4;C7",float),
                  ("sex_distrib","sex_dist","G4;G5",float),
                  ("smoking_status_distrib\tfemale",
                        ("ss",(len(GENDERS),len(SS), len(AGE_BRACKETS)), FEMALE),
                        "C12:Q14", float),
                  ("smoking_status_distrib\tmale",
                       ("ss",(len(GENDERS),len(SS), len(AGE_BRACKETS)), MALE),
                       "C18:Q20", float),
                  ("smoking_intensity_distrib\tfemale",
                       ("si",(len(GENDERS),len(SI), len(AGE_BRACKETS)), FEMALE),
                       "C26:Q28", float),
                  ("smoking_intensity_distrib\tmale",
                       ("si",(len(GENDERS),len(SI), len(AGE_BRACKETS)), MALE),
                       "C32:Q34", float),
                  ("former_smoker_quit_time_init_dist", "quit_dist","C38;C39",float),
                  ("event_prev","event_prev",add_prev,float)
                  
    ]

INPUTS['nathist'] = [("ns_lifetable", "ns_lifetable", "C6;D106",float),
                     ("xs_mort_usemult","xs_mort_usemult","I5",bool),
                     ("xs_mort_mult\tlight",
                          ("xs_mort_mult",(len(SI),len(GENDERS),len(AGE_BRACKETS)),SI_LIGHT),
                          "J15;K29",float),
                     ("xs_mort_mult\tmoderate",
                          ("xs_mort_mult",(len(SI),len(GENDERS),len(AGE_BRACKETS)),SI_MODERATE),
                          "N15;O29",float),
                     ("xs_mort_mult\theavy",
                          ("xs_mort_mult",(len(SI),len(GENDERS),len(AGE_BRACKETS)),SI_HEAVY),
                          "R15;S29",float),
                     ("xs_agequit_mult\tlight",
                          ("xs_agequit_mult",(len(SI),len(GENDERS),len(AGE_BRACKETS)),SI_LIGHT),
                          "J40;K54",float),
                     ("xs_agequit_mult\tmoderate",
                          ("xs_agequit_mult",(len(SI),len(GENDERS),len(AGE_BRACKETS)),SI_MODERATE),
                          "N40;O54",float),
                     ("xs_agequit_mult\theavy",
                          ("xs_agequit_mult",(len(SI),len(GENDERS),len(AGE_BRACKETS)),SI_HEAVY),
                          "R40;S54",float),
                     ("xs_lifetable\tlight",
                          ("xs_lifetable",(len(SI),len(GENDERS),101),SI_LIGHT),
                          "J61;K161",float),
                     ("xs_lifetable\tmoderate",
                          ("xs_lifetable",(len(SI),len(GENDERS),101),SI_MODERATE),
                          "N61;O161",float),
                     ("xs_lifetable\theavy",
                          ("xs_lifetable",(len(SI),len(GENDERS),101),SI_HEAVY),
                          "R61;S161",float),
                     ("nathist_transition","trans","W6:X6",int),
                     ("cs_mort_usemult","cs_mort_usemult","AG5",bool),
                     ("cs_mort_mult\tlight",
                          ("cs_mort_mult",(len(SI),len(GENDERS),len(AGE_BRACKETS)),SI_LIGHT),
                          "AH15;AI29",float),
                     ("cs_mort_mult\tmoderate",
                          ("cs_mort_mult",(len(SI),len(GENDERS),len(AGE_BRACKETS)),SI_MODERATE),
                          "AL15;AM29",float),
                     ("cs_mort_mult\theavy",
                          ("cs_mort_mult",(len(SI),len(GENDERS),len(AGE_BRACKETS)),SI_HEAVY),
                          "AP15;AQ29",float),
                     ("cs_lifetable\tlight",
                          ("cs_lifetable",(len(SI),len(GENDERS),101),SI_LIGHT),
                          "AH38;AI138",float),
                     ("cs_lifetable\tmoderate",
                          ("cs_lifetable",(len(SI),len(GENDERS),101),SI_MODERATE),
                          "AL38;AM138",float),
                     ("cs_lifetable\theavy",
                          ("cs_lifetable",(len(SI),len(GENDERS),101),SI_HEAVY),
                          "AP38;AQ138",float), 
    ]

INPUTS['smoking'] = [("prob_start_smoking","start_prob","C6;D106", float),
                    ("prob_quit_smoking","quit_prob","K6;L106", float),
                    ("event_quit_prob", "event_quit_prob", "P8;P12",float),
                    ("event_quit_duration", "event_quit_duration", "Q8;Q12",int),
                    ("comp_quit_prob", "comp_quit_prob", "P18;P22",float),
                    ("comp_quit_duration", "comp_quit_duration", "Q18;Q22",int),

                    ("smoking_relapse_coefficients\tlight",
                        ("relapse_coeffs", (len(SI), len(AGE_BRACKETS),2), SI_LIGHT),
                        "W7:X21", float),
                    ("smoking_relapse_coefficients\tmoderate",
                        ("relapse_coeffs", (len(SI), len(AGE_BRACKETS),2), SI_MODERATE),
                        "AA7:AB21", float),
                    ("smoking_relapse_coefficients\theavy",
                        ("relapse_coeffs", (len(SI), len(AGE_BRACKETS),2), SI_HEAVY),
                        "AE7:AF21", float),
    ]

INPUTS['pre_events'] = [("pevent_inc_baseline\tfemale",
                            ("pevent_inc_baseline",(len(GENDERS),NUM_EVENTS, len(AGE_BRACKETS)), FEMALE),
                            "C6:Q10", float),
                        ("pevent_inc_baseline\tmale",
                            ("pevent_inc_baseline",(len(GENDERS),NUM_EVENTS, len(AGE_BRACKETS)), MALE),
                            "C14:Q18", float),
                        ("pevent_inc_mult_curr_smokers", "pevent_inc_mult_curr", "C24:E28",float),
                        ("pevent_inc_mult_ex_smokers\tlight",
                            ("pevent_inc_mult_ex",(len(SI),NUM_EVENTS, len(AGE_BRACKETS)), SI_LIGHT),
                            "C35:Q39", float),
                        ("pevent_inc_mult_ex_smokers\tmoderate",
                            ("pevent_inc_mult_ex",(len(SI),NUM_EVENTS, len(AGE_BRACKETS)), SI_MODERATE),
                            "C43:Q47", float),
                        ("pevent_inc_mult_ex_smokers\theavy",
                            ("pevent_inc_mult_ex",(len(SI),NUM_EVENTS, len(AGE_BRACKETS)), SI_HEAVY),
                            "C51:Q55", float),
                        ("pevent_inc_transition", "pevent_inc_trans", "C61:D65", int),

                        ("pevent_mort_mult\tfemale",
                            ("pevent_mort_mult",(len(GENDERS),NUM_EVENTS, len(AGE_BRACKETS)), FEMALE),
                            "V6:AJ10", float),
                        ("pevent_mort_mult\tmale",
                            ("pevent_mort_mult",(len(GENDERS),NUM_EVENTS, len(AGE_BRACKETS)), MALE),
                            "V14:AJ18", float),
                        ("pevent_to_event_prob\tfemale",
                            ("pevent_to_event_prob",(len(GENDERS),NUM_EVENTS, len(AGE_BRACKETS)), FEMALE),
                            "V24:AJ28", float),
                        ("pevent_to_event_prob\tmale",
                            ("pevent_to_event_prob",(len(GENDERS),NUM_EVENTS, len(AGE_BRACKETS)), MALE),
                            "V32:AJ36", float),

                        ("screening_background_prob\tfemale",
                             ("screening_background_prob",(len(GENDERS),NUM_EVENTS,len(AGE_BRACKETS)),FEMALE),
                             "AO6:BC10", float),
                        ("screening_background_prob\tmale",
                             ("screening_background_prob",(len(GENDERS),NUM_EVENTS,len(AGE_BRACKETS)),MALE),
                             "AO14:BC18", float),
                        ("screening_sensitivity\tfemale",
                             ("screening_sensitivity",(len(GENDERS),NUM_EVENTS),FEMALE),
                             "AO24;AO28", float),
                        ("screening_sensitivity\tmale",
                             ("screening_sensitivity",(len(GENDERS),NUM_EVENTS),MALE),
                             "AO32;AO36", float),
                        ("screening_specificity\tfemale",
                             ("screening_specificity",(len(GENDERS),NUM_EVENTS),FEMALE),
                             "AO42;AO46", float),
                        ("screening_specificity\tmale",
                             ("screening_specificity",(len(GENDERS),NUM_EVENTS),MALE),
                             "AO50;AO54", float),
                        ("screening_regular_start_age", "screening_regular_start_age", "AP57:AT59", int),
                        ("screening_regular_interval", "screening_regular_interval", "AP60:AT60", int),
                        ("screening_regular_max", "screening_regular_max", "AP61:AT61", int),
                        ("screening_regular_prob_skip", "screening_regular_prob_skip", "AP62:AT62", float),
                        ("screening_conf_delay", "screening_conf_delay", "AP66:AT66", int),
                        ("screening_conf_mort", "screening_conf_mort", "AP67:AT67", float),
                        ("screening_outcome_neutral", "screening_outcome_neutral", "AP71:AT71", float),
                        ("screening_outcome_event_mult", "screening_outcome_event_mult", "AP72:AT72", float),
                        ("screening_outcome_pevent_mort_mult", "screening_outcome_pevent_mort_mult", "AP73:AT73", float),
                        ("screening_outcome_proph", "screening_outcome_proph", "AP74:AT74", int),
                        ("screening_outcome_intervention", "screening_outcome_intervention", "AP75:AT75", bool),
                        ]

INPUTS['events'] = [("event_inc_baseline\tfemale",
                        ("event_inc_baseline",(len(GENDERS),NUM_EVENTS, len(AGE_BRACKETS)), FEMALE),
                        "C6:Q10", float),
                    ("event_inc_baseline\tmale",
                        ("event_inc_baseline",(len(GENDERS),NUM_EVENTS, len(AGE_BRACKETS)), MALE),
                        "C14:Q18", float),
                    ("event_inc_mult_curr_smokers", "event_inc_mult_curr", "C24:E28",float),
                    ("event_inc_mult_ex_smokers\tlight",
                        ("event_inc_mult_ex",(len(SI),NUM_EVENTS, len(AGE_BRACKETS)), SI_LIGHT),
                        "C35:Q39", float),
                    ("event_inc_mult_ex_smokers\tmoderate",
                        ("event_inc_mult_ex",(len(SI),NUM_EVENTS, len(AGE_BRACKETS)), SI_MODERATE),
                        "C43:Q47", float),
                    ("event_inc_mult_ex_smokers\theavy",
                        ("event_inc_mult_ex",(len(SI),NUM_EVENTS, len(AGE_BRACKETS)), SI_HEAVY),
                        "C51:Q55", float),
                    ("event_inc_transition", "event_inc_trans", "C61:D65", int),
                    ("event_prob_death","event_prob_death","C70;C74",float),

                    ("event_comp_baseline\tfemale",
                        ("event_comp_baseline",(len(GENDERS),NUM_EVENTS, len(AGE_BRACKETS)), FEMALE),
                        "V6:AJ10", float),
                    ("event_comp_baseline\tmale",
                        ("event_comp_baseline",(len(GENDERS),NUM_EVENTS, len(AGE_BRACKETS)), MALE),
                        "V14:AJ18", float),
                    ("event_comp_mult_curr_smokers", "event_comp_mult_curr", "V24:X28",float),
                    ("event_comp_mult_ex_smokers\tlight",
                        ("event_comp_mult_ex",(len(SI),NUM_EVENTS, len(AGE_BRACKETS)), SI_LIGHT),
                        "V35:AJ39", float),
                    ("event_comp_mult_ex_smokers\tmoderate",
                        ("event_comp_mult_ex",(len(SI),NUM_EVENTS, len(AGE_BRACKETS)), SI_MODERATE),
                        "V43:AJ47", float),
                    ("event_comp_mult_ex_smokers\theavy",
                        ("event_comp_mult_ex",(len(SI),NUM_EVENTS, len(AGE_BRACKETS)), SI_HEAVY),
                        "V51:AJ55", float),
                    ("event_comp_transition", "event_comp_trans", "V61:W65", int),
                    ("event_comp_prob_death","event_comp_prob_death","V70;V74",float),
                    ]


INPUTS['prophs'] = [("proph_names","proph_names","C6;C15", str),
                    ("proph_eff_events", "eff_events", "D6:H15", float),
                    ("proph_eff_comp", "eff_comp", "D21:H30", float),
                    ("proph_tox_prob", "tox_prob", "D36;D45", float),
                    ("proph_tox_dth_prob", "tox_dth_prob", "E36;E45", float),
                    ("proph_stop_on_tox", "stop_on_tox", "F36;F45", bool),
                    ("proph_allow_restart_on_tox", "allow_restart_on_tox", "G36;G45", bool),

                    ("proph_start_prob_month\tlight",
                         ("start_prob_month", (len(SI),len(GENDERS),NUM_PROPHS), SI_LIGHT),
                         "N7;O16", float),
                    ("proph_start_prob_month\tmoderate",
                         ("start_prob_month", (len(SI),len(GENDERS),NUM_PROPHS), SI_MODERATE),
                         "P7;Q16", float),
                    ("proph_start_prob_month\theavy",
                         ("start_prob_month", (len(SI),len(GENDERS),NUM_PROPHS), SI_HEAVY),
                         "R7;S16", float),
                    ("proph_start_prob_init\tlight",
                         ("start_prob_init", (len(SI),len(GENDERS),NUM_PROPHS), SI_LIGHT),
                         "N23;O32", float),
                    ("proph_start_prob_init\tmoderate",
                         ("start_prob_init", (len(SI),len(GENDERS),NUM_PROPHS), SI_MODERATE),
                         "P23;Q32", float),
                    ("proph_start_prob_init\theavy",
                         ("start_prob_init", (len(SI),len(GENDERS),NUM_PROPHS), SI_HEAVY),
                         "R23;S32", float),
                    ("proph_start_age_mult\tfemale",
                         ("start_age_mult", (len(GENDERS),NUM_PROPHS,len(AGE_BRACKETS)), FEMALE),
                         "N38:AB47", float),
                    ("proph_start_age_mult\tmale",
                         ("start_age_mult", (len(GENDERS),NUM_PROPHS,len(AGE_BRACKETS)), MALE),
                         "N51:AB60", float),
                    ("proph_start_hist_mult\tfemale",
                         ("event_hist_mult", (len(GENDERS),NUM_PROPHS,NUM_EVENTS), FEMALE),
                         "N66:R75", float),
                    ("proph_start_hist_mult\tmale",
                         ("event_hist_mult", (len(GENDERS),NUM_PROPHS,NUM_EVENTS), MALE),
                         "N79:R88", float),
                    ("proph_start_ss_mult\tfemale",
                         ("ss_mult", (len(GENDERS),NUM_PROPHS,len(SS)), FEMALE),
                         "N94:P103", float),
                    ("proph_start_ss_mult\tmale",
                         ("ss_mult", (len(GENDERS),NUM_PROPHS,len(SS)), MALE),
                         "N107:P116", float),

                    ("proph_stop_prob_month\tlight",
                         ("stop_prob_month", (len(SI),len(GENDERS),NUM_PROPHS), SI_LIGHT),
                         "AH7;AI16", float),
                    ("proph_stop_prob_month\tmoderate",
                         ("stop_prob_month", (len(SI),len(GENDERS),NUM_PROPHS), SI_MODERATE),
                         "AJ7;AK16", float),
                    ("proph_stop_prob_month\theavy",
                         ("stop_prob_month", (len(SI),len(GENDERS),NUM_PROPHS), SI_HEAVY),
                         "AL7;AM16", float),
                    ]



INPUTS['intervention'] = [("int_names","int_names","C6;C15", str),
                    ("int_quit_mult", "quit_mult", "D6;D15", float),
                    ("int_relapse_mult", "relapse_mult", "E6;E15", float),
                    ("int_tox_prob", "tox_prob", "D21;D30", float),
                    ("int_tox_dth_prob", "tox_dth_prob", "E21;E30", float),
                    ("int_stop_on_tox", "stop_on_tox", "F21;F30", bool),
                    ("int_allow_restart_on_tox", "allow_restart_on_tox", "G21;G30", bool),

                    ("int_start_prob_month\tlight",
                         ("start_prob_month", (len(SI),len(GENDERS),NUM_PROPHS), SI_LIGHT),
                         "N7;O16", float),
                    ("int_start_prob_month\tmoderate",
                         ("start_prob_month", (len(SI),len(GENDERS),NUM_PROPHS), SI_MODERATE),
                         "P7;Q16", float),
                    ("int_start_prob_month\theavy",
                         ("start_prob_month", (len(SI),len(GENDERS),NUM_PROPHS), SI_HEAVY),
                         "R7;S16", float),
                    ("int_start_prob_init\tlight",
                         ("start_prob_init", (len(SI),len(GENDERS),NUM_PROPHS), SI_LIGHT),
                         "N23;O32", float),
                    ("int_start_prob_init\tmoderate",
                         ("start_prob_init", (len(SI),len(GENDERS),NUM_PROPHS), SI_MODERATE),
                         "P23;Q32", float),
                    ("int_start_prob_init\theavy",
                         ("start_prob_init", (len(SI),len(GENDERS),NUM_PROPHS), SI_HEAVY),
                         "R23;S32", float),
                    ("int_start_age_mult\tfemale",
                         ("start_age_mult", (len(GENDERS),NUM_PROPHS,len(AGE_BRACKETS)), FEMALE),
                         "N38:AB47", float),
                    ("int_start_age_mult\tmale",
                         ("start_age_mult", (len(GENDERS),NUM_PROPHS,len(AGE_BRACKETS)), MALE),
                         "N51:AB60", float),
                    ("int_start_hist_mult\tfemale",
                         ("event_hist_mult", (len(GENDERS),NUM_PROPHS,NUM_EVENTS), FEMALE),
                         "N66:R75", float),
                    ("int_start_hist_mult\tmale",
                         ("event_hist_mult", (len(GENDERS),NUM_PROPHS,NUM_EVENTS), MALE),
                         "N79:R88", float),

                    ("int_stop_prob_month\tlight",
                         ("stop_prob_month", (len(SI),len(GENDERS),NUM_PROPHS), SI_LIGHT),
                         "AH7;AI16", float),
                    ("int_stop_prob_month\tmoderate",
                         ("stop_prob_month", (len(SI),len(GENDERS),NUM_PROPHS), SI_MODERATE),
                         "AJ7;AK16", float),
                    ("int_stop_prob_month\theavy",
                         ("stop_prob_month", (len(SI),len(GENDERS),NUM_PROPHS), SI_HEAVY),
                         "AL7;AM16", float),
                    ("int_stop_on_quit", "stop_on_quit", "AH22;AH31", bool),
                    ("int_stop_abst_duration", "stop_abst_duration", "AI22;AI31", int),
                    ("int_duration", "duration", "AJ22;AJ31", int),
                    ]


INPUTS['costs'] = [("costs_event_init","event_init","C5;C9", float),
                   ("costs_event_month","event_month","D5;D9", float),
                   ("costs_comp","comp","C14;C18", float),
                   ("costs_proph_init","proph_init","C23;C32", float),
                   ("costs_proph_month","proph_month","D23;D32", float),
                   ("costs_int_init","int_init","C37;C46", float),
                   ("costs_int_month","int_month","D37;D46", float),

                   ("costs_bkgd_trans", "bkgd_trans","K5",int),
                   ("costs_bkgd\tfemale",
                         ("bkgd", (len(GENDERS),len(SS),len(AGE_BRACKETS)), FEMALE),
                         "I8:W10", float),
                   ("costs_bkgd\tmale",
                         ("bkgd", (len(GENDERS),len(SS),len(AGE_BRACKETS)), MALE),
                         "I14:W16", float),
                   ("costs_screen","screen","I22;I26", float),
                   ("costs_screen_pos","screen_pos","J22;J26", float),
                   ("costs_screen_neg","screen_neg","K22;K26", float),
                   ("costs_screen_conf","screen_conf","L22;L26", float),
                   ("costs_screen_detected","screen_detected","M22;M26", float),
                    ]

INPUTS['qol'] = [("enable_qol","enable_qol","B5", bool),
                 ("qol_base","base","C11:E13", float),
                 ("qol_quit","quit","C19;C21", float),
                 ("qol_quit_duration","quit_duration","D19;D21", int),

                 ("qol_screen","screen","J6;J10", float),
                 ("qol_screen_wait_conf","screen_wait_conf","K6;K10", float),
                 ("qol_screen_conf","screen_conf","L6;L10", float),
                 ("qol_screen_undet","screen_undet","M6;M10", float),
                 ("qol_screen_det","screen_det","N6;N10", float),
                 ("qol_event_init","event_init","J16;J20", float),
                 ("qol_event_month","event_month","K16;K20", float),
                 ("qol_comp","comp","J26;J30", float),

                 ("qol_proph_tox","proph_tox","S6;S15", float),
                 ("qol_int_tox","int_tox","S21;S30", float),
                 ]
#gets a value from a range of cells
#can be a single cell (E5) or range (E5:F6)
#convention is to use : for rowslices or ; for colslices
def xl_get_range(sheet, address, tp):
    #single cell
    if ":" not in address and ";" not in address:
        return tp(sheet.cell(*xl_cell_to_rowcol(address)).value)
    #row slice
    if ":" in address:
        start, end = address.split(":")
        startrow, startcol = xl_cell_to_rowcol(start)
        endrow, endcol = xl_cell_to_rowcol(end)
        if startrow == endrow:
            return np.array([tp(cell) for cell in sheet.row_values(startrow, startcol, endcol+1)])
        else:
            values = []
            for row in range(startrow, endrow+1):
                values.append([tp(cell) for cell in sheet.row_values(row, startcol, endcol+1)])

            return np.array(values)
    #col slice
    if ";" in address:
        start, end = address.split(";")
        startrow, startcol = xl_cell_to_rowcol(start)
        endrow, endcol = xl_cell_to_rowcol(end)
        if startcol == endcol:
            return np.array([tp(cell) for cell in sheet.col_values(startcol, startrow, endrow+1)])
        else:
            values = []
            for col in range(startcol, endcol+1):
                values.append([tp(cell) for cell in sheet.col_values(col, startrow, endrow+1)])
                
            return np.array(values)

def text_get_range(textfile, address, label, tp):

    if tp == bool:
        tp = lambda x:x=="True"
    # single cell
    if ":" not in address and ";" not in address:
        val = re.search(label+"\t(.*)\n",textfile).group(1)
        return tp(val)

    # row slice
    if ":" in address:
        start, end = address.split(":")
        startrow, startcol = xl_cell_to_rowcol(start)
        endrow, endcol = xl_cell_to_rowcol(end)
        if startrow == endrow:
            val = re.search(label+"\t(.*)\n",textfile).group(1)
            return np.array([tp(_) for _ in val.split("\t")])
        else:
            values = []
            pattern = label+"\t"
            for row in range(startrow, endrow + 1):
                pattern+="(.*)\n"

            val = re.search(pattern,textfile)
            for row in val.groups():
                values.append([tp(_) for _ in row.split("\t")])
            return np.array(values)
    # col slice
    if ";" in address:
        start, end = address.split(";")
        startrow, startcol = xl_cell_to_rowcol(start)
        endrow, endcol = xl_cell_to_rowcol(end)
        if startcol == endcol:
            val = re.search(label+"\t(.*)\n",textfile).group(1)
            return np.array([tp(_) for _ in val.split("\t")])
        else:
            values = []
            pattern = label+"\t"
            for col in range(startcol, endcol + 1):
                pattern+="(.*)\n"

            val = re.search(pattern,textfile)

            for col in val.groups():
                values.append([tp(_) for _ in col.split("\t")])
            return np.array(values)

def text_write_range(textfile, address, label, val):
    #single cell
    if ":" not in address and ";" not in address:
        textfile.write("{}\t{}".format(label, val))
    #row or col slice
    else:
        textfile.write("{}\t".format(label))
        if len(val.shape)==1:
            textfile.write("\t".join([str(_) for _ in val]))
        else:
            for i in range(val.shape[0]):
                if i!=0:
                    textfile.write("\n")
                textfile.write("\t".join([str(_) for _ in val[i]]))

    textfile.write("\n")

#converts excel address to row,col
def xl_cell_to_rowcol(cell):
    collet,rownum = re.match("([a-zA-Z]*)(\d*)",cell).groups()
    col = 0
    for let in collet.upper():
        col*=26
        col+=ord(let)-64
    col-=1
    return int(rownum)-1,col
        
        
#class for tab in excel sheet
class Tab(object):
    def __init__(self, sheetname):
        self.sheetname = sheetname
        self.inputnames = []
    def add_input(self, varname, value):
        if isinstance(varname, tuple):
            base, shape, indices = varname
            
            try:
                getattr(self,base)
            except AttributeError:
                setattr(self,base, np.ndarray(shape))
                self.inputnames.append(base)

            getattr(self,base)[indices] = value
        else:

            setattr(self, varname, value)
            self.inputnames.append(varname)
    def __repr__(self):
        return str(self.inputnames)
#class for inputs.  Each tab has it's own var.
#inputs can be accessed using . notation e.g. inputs.sim.runsize
class Inputs(object):
    def __init__(self):
        #create object for each tab
        self.tabnames = []

    def load_excel(self, filepath):
        self.delete_inputs()
        for sheetname, varname  in TABS:
            setattr(self, varname, Tab(sheetname))
            self.tabnames.append(varname)
        wb = xlrd.open_workbook(filepath)
        for sheetname, sheetvar in TABS:
            inputtab = INPUTS[sheetvar]
            sheet = wb.sheet_by_name(sheetname)
            for label, varname, address, tp in inputtab:
                if isinstance(address, collections.Callable):
                    #its a function
                    value = address(sheet)
                else:
                    value = xl_get_range(sheet, address, tp)
                getattr(self,sheetvar).add_input(varname,value)

    #saves inputs to text file
    def save_txt(self, filepath):
        with open(filepath,'w') as fwrite:
            for sheetname, sheetvar in TABS:
                inputtab = INPUTS[sheetvar]
                for label, varname, address, tp in inputtab:
                    if isinstance(address, collections.Callable):
                        fname = address.__name__
                        fname+="_write_text"
                        globals()[fname](fwrite,self)
                    else:
                        if isinstance(varname, tuple):
                            base, shape, indices = varname
                            value = getattr(getattr(self, sheetvar), base)[indices]
                        else:
                            value = getattr(getattr(self,sheetvar),varname)
                        text_write_range(fwrite,address,label,value)

    #loads inputs from text file
    def load_txt(self, filepath):
        self.delete_inputs()
        for sheetname, varname  in TABS:
            setattr(self, varname, Tab(sheetname))
            self.tabnames.append(varname)

        with open(filepath) as fread:
            textfile = fread.read()

        for sheetname, sheetvar in TABS:
            inputtab = INPUTS[sheetvar]
            for label, varname, address, tp in inputtab:
                if isinstance(address, collections.Callable):
                    #its a function
                    fname = address.__name__
                    fname += "_load_text"
                    value =globals()[fname](textfile, self)
                else:
                    value = text_get_range(textfile, address,label, tp)
                getattr(self,sheetvar).add_input(varname,value)

    def delete_inputs(self):
        self.tabnames = []
        for sheetname, varname in TABS:
            if hasattr(self, varname):
                delattr(self, varname)
    def __repr__(self):
        return str(self.tabnames)
if __name__ == "__main__":
    excelfile = "input file/Smoking Model Inputs.xlsm"
    i = Inputs()
    i.load_txt("test.smin")
    #i.load_excel(excelfile)
    #i.save_txt("test.smin")

        
