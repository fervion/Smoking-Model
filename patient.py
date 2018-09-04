"""
The patient object for smoking model

@author: Taige Hou(thou1@partners.org)
"""

from enums import *
from inputs import *
from outputs import *
from math import exp

######################################################################
class Patient(object):
    def __init__(self, pid,inputs, outputs, disc_outputs, monthly_outputs, traces = None):
        self.inputs = inputs
        self.outputs = outputs
        self.disc_outputs = disc_outputs
        self.monthly_outputs = monthly_outputs
        #current monthly output object
        self.curr_mo= None
        self.pid = pid
        self.traces = traces
        self.trace_text = ""
        #calendar month of run
        self.month = 0
        self.isalive = True
        self.causeofdeath = None

        self.update_qol = self.inputs.qol.enable_qol

        #month of next regular screening
        self.month_of_reg_screen = [None for i in range(NUM_EVENTS)]
        self.num_reg_screens = [0 for i in range(NUM_EVENTS)]
        self.month_of_conf = [None for i in range(NUM_EVENTS)]
        self.is_detected_pevent = [False for i in range(NUM_EVENTS)]
        
        #mortality risks
        self.mort_risks = []

        #smoking
        self.ever_quit = False
        self.ever_relapse = False
        
        #interventions
        self.has_int = [False for i in range(NUM_INTERVENTIONS)]
        self.month_start_int = [None for i in range(NUM_INTERVENTIONS)]
        self.has_int_tox_hist = [False for i in range(NUM_INTERVENTIONS)]
        self.ever_start_any_int = False
        
        #prophs
        self.has_proph = [False for i in range(NUM_PROPHS)]
        self.has_proph_tox_hist = [False for i in range(NUM_PROPHS)]

        #events
        self.events = [EVT_NONE for i in range(NUM_EVENTS)]
        self.month_of_event = [None for i in range(NUM_EVENTS)]
        self.month_of_comp = [None for i in range(NUM_EVENTS)]
        
        self.init_patient()
    #initialize patient
    def init_patient(self):
        init = self.inputs.init

        self.event_names = self.inputs.sim.event_names
        self.int_names = self.inputs.intervention.int_names
        self.proph_names = self.inputs.prophs.proph_names
        self.death_names=["Nat Hist","Old Age", "Proph Tox", "Intervention Tox", "Confirmatory Test",]
        self.death_names.extend(["{0}".format(self.event_names[event]) for event in range(NUM_EVENTS)])
        self.death_names.extend(["{0} Complication".format(self.event_names[event]) for event in range(NUM_EVENTS)])

        self.disc_mult_month = pow(self.inputs.sim.disc_rate_year+1,1/12.0)
        self.disc_factor = 1.0
        self.gender = draw_dist(init.sex_dist)
        self.age = int(round(draw_trunc_norm(*init.age_dist)))
        self.agecat = get_age_cat(AGE_BRACKETS, self.age)
        self.ss = draw_dist(init.ss[self.gender,:,self.agecat])
        self.si = draw_dist(init.si[self.gender,:,self.agecat])

        #roll for time since quit
        self.month_of_quit = None
        if self.ss == SS_FORMER:
            self.month_of_quit = -1*int(round(draw_trunc_norm(*(init.quit_dist),low=0,high=self.age)))
            self.agequit = self.age - (self.month - self.month_of_quit)
            self.agequit_cat = get_age_cat(AGE_BRACKETS, self.agequit)

        
        #draw event prev
        for event in range(NUM_EVENTS):
            prob = init.event_prev[event,self.gender,self.ss,self.si,self.agecat]
            if random() <= prob:
                self.get_event(event, True)

        #roll for start proph at init
        prophs = self.inputs.prophs
        for i in range(NUM_PROPHS):
            prob_proph = prophs.start_prob_init[self.si,self.gender,i]
            prob_proph *= prophs.start_age_mult[self.gender, i, self.agecat]
            evts = [j for j in range(NUM_EVENTS) if self.events[j] == EVT_FULL]
            prob_proph *= prophs.event_hist_mult[self.gender, i, evts].prod()
            prob_proph *= prophs.ss_mult[self.gender, i, self.ss]
                
            if random() <= prob_proph:
                self.start_proph(i, True)

        #roll for start intervention at init
        intv = self.inputs.intervention
        for i in range(NUM_INTERVENTIONS):
            prob_int = intv.start_prob_init[self.si,self.gender,i]
            prob_int *= intv.start_age_mult[self.gender, i, self.agecat]
            evts = [j for j in range(NUM_EVENTS) if self.events[j] == EVT_FULL]
            mults = intv.event_hist_mult[self.gender, i, evts]
            if mults.size:
                prob_int *= mults.max()
            
            if random() <= prob_int:
                self.start_intervention(i, True)

        #set next regular screening month
        for event in range(NUM_EVENTS):
            if self.inputs.pre_events.screening_regular_max[event] > 0:
                start_age = self.inputs.pre_events.screening_regular_start_age[self.ss][event]*12
                if start_age >= 0:
                    self.month_of_reg_screen[event] = start_age - self.age

        #add outputs
        self.add_out("init_dist_smoking",(self.ss, self.si),1)
        self.add_out("init_dist_gender",(self.gender,),1)
        self.add_out("init_age",None, self.age)

        #prints to trace file if it exists
        if self.traces is not None:
            self.trace("\n\nBEGIN PATIENT {0}".format(self.pid), False)
            self.trace("\n\tgender: {0}, init age: {1} mths".format(GENDER_STRS[self.gender], self.age), False)
            self.trace("\n\tsmoking status: {0}, smoking intensity: {1}".format(SS_STRS[self.ss],SI_STRS[self.si]), False)
            if self.ss == SS_FORMER:
                self.trace("\n\tmonth of quit: {0}".format(self.month_of_quit))
            self.trace("\n\tevents: {0}".format(",".join(
                [self.event_names[event] for event in range(NUM_EVENTS) if self.events[event]==EVT_FULL])), False)
            self.trace("\n\tinterventions: {0}".format(",".join(
                [self.int_names[intv] for intv in range(NUM_INTERVENTIONS) if self.has_int[intv]])), False)
            self.trace("\n\tprophs: {0}".format(",".join(
                [self.proph_names[proph] for proph in range(NUM_PROPHS) if self.has_proph[proph]])), False)
    def update_start_month(self):
        #clear mortality risks
        self.mort_risks = []
        #calculate age_category
        self.agecat = get_age_cat(AGE_BRACKETS, self.age)
        self.qol = 1.0

        #add month for output
        try:
            self.curr_mo = self.monthly_outputs[self.month]
        except IndexError:
            self.curr_mo = self.monthly_outputs.add_month()
        
    def update_interventions(self):
        if self.ss == SS_NEVER:
            return
        
        intv = self.inputs.intervention
        
        for i in range(NUM_INTERVENTIONS):
            #those who have intervention
            if self.has_int[i]:
                
                #Check for int stop
                stop_int = False
                
                #abstinance condition
                if self.ss==SS_FORMER:
                    quit_duration = self.month - self.month_of_quit
                    abst_cond = intv.stop_abst_duration[i]
                    if abst_cond>=0 and quit_duration >= abst_cond:
                        stop_int = True
                        
                dur_cond = intv.duration[i]
                int_duration = self.month - self.month_start_int[i]
                if not stop_int and dur_cond >=0 and int_duration >= dur_cond:
                    stop_int = True

                if not stop_int:
                    prob_stop = intv.stop_prob_month[self.si, self.gender, i]        
                    if random() <= prob_stop:
                        stop_int = True

                if stop_int:
                    self.stop_intervention(i)
            #check for int start
            if self.ss==SS_CURRENT and not self.has_int[i]:
                if not self.has_int_tox_hist[i] or intv.allow_restart_on_tox[i]:
                    prob_start = intv.start_prob_month[self.si, self.gender, i]
                    prob_start *= intv.start_age_mult[self.gender, i, self.agecat]
                    evts = [j for j in range(NUM_EVENTS) if self.events[j] == EVT_FULL]
                    mults = intv.event_hist_mult[self.gender, i, evts]
                    if mults.size:
                        prob_start *= mults.max()                    
                    if random() <= prob_start:
                        self.start_intervention(i)

            #toxicity
            if self.has_int[i] and random() <= intv.tox_prob[i]:
                self.set_int_tox(i)
                #add tox death prob
                self.add_mort_risk(DTH_TOX_INT,intv.tox_dth_prob[i])
                #stop intervention if set
                if intv.stop_on_tox[i]:
                    self.stop_intervention(i)
                
    def update_smoke_start(self):
        #updater for starting smoking
        if self.ss != SS_NEVER:
            return
        ageyrs = self.age // 12
        probstart = self.inputs.smoking.start_prob[self.gender,ageyrs]
        #roll for start
        if random() <= probstart:
            self.start_smoking()

    def update_smoke_quit(self):
        #updater for stoping smoking
        if self.ss != SS_CURRENT:
            return

        smoking = self.inputs.smoking
        ageyrs = self.age // 12


        evts = [i for i in range(NUM_EVENTS) if self.events[i] == EVT_FULL and
                self.month - self.month_of_event[i]<=smoking.event_quit_duration[i]]
        comps = [i for i in range(NUM_EVENTS) if self.month_of_comp[i] is not None and
                self.month - self.month_of_comp[i]<=smoking.comp_quit_duration[i]]

        probs = np.concatenate((smoking.event_quit_prob[evts],smoking.comp_quit_prob[comps]))
        if probs.size:
            probquit= max(probs)
        else:
            probquit = smoking.quit_prob[self.gender, ageyrs]

        int_mults = self.inputs.intervention.quit_mult[self.has_int]


        if int_mults.size:
            probquit *= int_mults.max()

        if random() <= probquit:
            self.quit_smoking()
    
    def update_smoke_relapse(self):
        #updater for smoking relapse
        if self.ss != SS_FORMER:
            return
        smoking = self.inputs.smoking

        time_since_quit = self.month - self.month_of_quit
        c,b = smoking.relapse_coeffs[self.si, self.agecat]
        if time_since_quit==0:
            probrelapse = 0
        else:
            probrelapse = c*exp(b*time_since_quit)

        int_mults = self.inputs.intervention.relapse_mult[self.has_int]

        if int_mults.size:
            probrelapse *= int_mults.max()

        if random() <= probrelapse:
            self.relapse_smoking()
            
    def update_prophs(self):
        prophs = self.inputs.prophs
        
        for i in range(NUM_PROPHS):
            #those with proph
            if self.has_proph[i]:
                
                #Check for proph stop
                prob_stop = prophs.stop_prob_month[self.si, self.gender, i]        
                if random() <= prob_stop:
                    self.stop_proph(i)

            #check for proph start
            if not self.has_proph[i]:
                if not self.has_proph_tox_hist[i] or prophs.allow_restart_on_tox[i]:
                    prob_start = prophs.start_prob_month[self.si, self.gender, i]
                    prob_start *= prophs.start_age_mult[self.gender, i, self.agecat]
                    evts = [j for j in range(NUM_EVENTS) if self.events[j] == EVT_FULL]
                    mults = prophs.event_hist_mult[self.gender, i, evts]
                    if mults.size:
                        prob_start *= mults.max()
                    prob_start *= prophs.ss_mult[self.gender, i, self.ss]
                    if random() <= prob_start:
                        self.start_proph(i)

            #toxicity
            if self.has_proph[i] and random() <= prophs.tox_prob[i]:
                self.set_proph_tox(i)
                #add tox death prob
                self.add_mort_risk(DTH_TOX_PROPH,prophs.tox_dth_prob[i])
                #stop proph if set
                if prophs.stop_on_tox[i]:
                    self.stop_proph(i)
                
    def update_pre_events(self):
        #roll for incidence of pre-events
        pevt = self.inputs.pre_events
        agecat = self.agecat
        
        for event in range(NUM_EVENTS):
            if self.events[event] != EVT_NONE:
                continue
            
            prob_event = pevt.pevent_inc_baseline[self.gender,event,agecat]
            cs_mult = pevt.pevent_inc_mult_curr[event,self.si]
            
            if self.ss == SS_CURRENT:
                prob_event *= cs_mult
            elif self.ss == SS_FORMER:
                trans = pevt.pevent_inc_trans[event]
                time_since_quit = self.month - self.month_of_quit
                xs_mult = pevt.pevent_inc_mult_ex[self.si,event,self.agequit_cat]
                #ex smokers
                if time_since_quit >= trans[1]:
                    prob_event *= xs_mult
                #current smokers
                elif time_since_quit <= trans[0]:
                    prob_event *= cs_mult
                #stopped smokers                    
                else:
                    prob_event *= np.interp(time_since_quit, trans, [cs_mult, xs_mult])
            #roll for pre-event
            if random() <= prob_event:
                self.get_pevent(event)

        for event in range(NUM_EVENTS):
            if self.events[event] == EVT_PRE:
                self.add_out("pevent_mths",event,1)
    def update_pevent_screening(self):
        #handles screening for pre events
        pevt = self.inputs.pre_events
        
        for event in range(NUM_EVENTS):
            if self.events[event] == EVT_FULL:
                continue            
            
            if self.month_of_conf[event] is None:
                #check if regular screening should be done
                has_reg_screen = False
                month_of_screen = self.month_of_reg_screen[event]

                if month_of_screen is not None and self.month >= month_of_screen:
                    #roll for prob attend
                    if random() > pevt.screening_regular_prob_skip[event]:
                        self.screen(SCREEN_REG, event)
                        self.num_reg_screens[event]+=1
                        has_reg_screen = True
                        
                    #schedule next reg screen if possible
                    if self.num_reg_screens[event] < pevt.screening_regular_max[event]:
                        next_screen = self.month+pevt.screening_regular_interval[event]
                    else:
                        next_screen = None
                    self.month_of_reg_screen[event] = next_screen

                #roll for background screening
                if not has_reg_screen and random() <= pevt.screening_background_prob[self.gender,event,self.agecat]:
                    self.screen(SCREEN_BACK, event)


            #awaiting conf test
            if self.month_of_conf[event] is not None:
                #add qol
                if self.update_qol:
                    self.add_qol(self.inputs.qol.screen_wait_conf[event])
                    
                #check for conf test
                if self.month==self.month_of_conf[event]:
                    self.conf_test(event)

            if self.events[event] == EVT_PRE:
                #add cost for being detected
                if self.is_detected_pevent[event]:
                    self.add_cost("cost_screening",event, self.inputs.costs.screen_detected[event], True)
                    #add qol
                    if self.update_qol:
                        self.add_qol(self.inputs.qol.screen_det[event])
                else:
                    if self.update_qol:
                        self.add_qol(self.inputs.qol.screen_undet[event])
                    
    def update_events(self):
        #updater for events
        inp = self.inputs.events
        agecat = self.agecat
        prophs = [i for i in range(NUM_PROPHS) if self.has_proph[i]]

        #for those who have event
        for event in range(NUM_EVENTS):
            if self.events[event] != EVT_FULL:
                continue
            
            #add monthly cost of event:
            self.add_cost("cost_event",event, self.inputs.costs.event_month[event], True)

            #add qol
            if self.update_qol:
                self.add_qol(self.inputs.qol.event_month[event])
                
            #roll for event complications
            prob_comp = inp.event_comp_baseline[self.gender, event, agecat]
            cs_mult = inp.event_comp_mult_curr[event, self.si]
            
            if self.ss == SS_CURRENT:
                prob_comp *= cs_mult
            elif self.ss == SS_FORMER:
                trans = inp.event_comp_trans[event]
                time_since_quit = self.month - self.month_of_quit
                xs_mult = inp.event_comp_mult_ex[self.si, event, self.agequit_cat]
                #ex smokers
                if time_since_quit >= trans[1]:
                    prob_comp *= xs_mult
                #current smokers
                elif time_since_quit <= trans[0]:
                    prob_comp *= cs_mult
                #stopped smokers
                else:
                    prob_comp *= np.interp(time_since_quit, trans, [cs_mult, xs_mult])
                    
            #modify by proph
            prob_comp *= self.inputs.prophs.eff_comp[prophs,event].prod()
            
            #roll for complication
            if random() <= prob_comp:
                self.get_event_comp(event)
                #roll for complication death
                prob_death = inp.event_comp_prob_death[event]
                self.add_mort_risk(globals()["DTH_EVENT_COMP_{0}".format(event)],prob_death)
                
        #roll for transitions from pre-event to event
        for event in range(NUM_EVENTS):
            if self.events[event] == EVT_PRE:
                prob_event = self.inputs.pre_events.pevent_to_event_prob[self.gender, event, agecat]
                #modify by detected status
                if self.is_detected_pevent[event]:
                    prob_event*=self.inputs.pre_events.screening_outcome_event_mult[event]
                if random() <= prob_event:
                    self.get_event(event)
                    
        #roll for incidence of event
        for event in range(NUM_EVENTS):
            if self.events[event] != EVT_NONE:
                continue
            
            prob_event = inp.event_inc_baseline[self.gender,event,agecat]
            cs_mult = inp.event_inc_mult_curr[event,self.si]
            
            if self.ss == SS_CURRENT:
                prob_event *= cs_mult
            elif self.ss == SS_FORMER:
                trans = inp.event_inc_trans[event]
                time_since_quit = self.month - self.month_of_quit
                xs_mult = inp.event_inc_mult_ex[self.si,event,self.agequit_cat]
                #ex smokers
                if time_since_quit >= trans[1]:
                    prob_event *= xs_mult
                #current smokers
                elif time_since_quit <= trans[0]:
                    prob_event *= cs_mult
                #stopped smokers                    
                else:
                    prob_event *= np.interp(time_since_quit, trans, [cs_mult, xs_mult])

            #modify by proph
            prob_event *= self.inputs.prophs.eff_events[prophs,event].prod()

            #roll for event
            if random() <= prob_event:
                self.get_event(event)
                #roll for event death
                prob_death = inp.event_prob_death[event]
                self.add_mort_risk(globals()["DTH_EVENT_{0}".format(event)],prob_death)

    def update_nathist(self):
        nh = self.inputs.nathist
        ageyrs = self.age//12
        agecat = self.agecat
        
        #accumulate pre-event multiplier
        pevent_mult = 1
        for event in range(NUM_EVENTS):
            if self.events[event] == EVT_PRE:
                pevent_mult *= self.inputs.pre_events.pevent_mort_mult[self.gender,event, agecat]
                if self.is_detected_pevent[event]:
                    pevent_mult *= self.inputs.pre_events.screening_outcome_pevent_mort_mult[event]

        ns_mort = nh.ns_lifetable[self.gender,ageyrs]
        #never smokers
        if (self.ss == SS_NEVER):
            self.add_mort_risk(DTH_NAT_HIST,ns_mort,pevent_mult)
            return
        
        #current smokers
        if nh.cs_mort_usemult:
            cs_mort = ns_mort*nh.cs_mort_mult[self.si,self.gender,agecat]
        else:
            cs_mort = nh.cs_lifetable[self.si,self.gender,ageyrs]

        if (self.ss == SS_CURRENT):
            self.add_mort_risk(DTH_NAT_HIST,cs_mort, pevent_mult)
            return
                
        #former smokers
        time_since_quit = self.month - self.month_of_quit
        trans = nh.trans
        #ex smokers
        if nh.xs_mort_usemult:
            xs_mort = ns_mort * nh.xs_mort_mult[self.si,self.gender,agecat]
        else:
            age_at_quit = self.age - time_since_quit
            xs_mort = nh.xs_agequit_mult[self.si,self.gender,agecat]*nh.xs_lifetable[self.si,self.gender,ageyrs]

        #ex smokers
        if time_since_quit >= trans[1]:
            self.add_mort_risk(DTH_NAT_HIST, xs_mort, pevent_mult)
        #current smokers
        elif time_since_quit <= trans[0]:
            self.add_mort_risk(DTH_NAT_HIST, cs_mort, pevent_mult)
        #stopped smokers
        else:
            self.add_mort_risk(DTH_NAT_HIST,np.interp(time_since_quit, trans, [cs_mort, xs_mort]), pevent_mult)

    def update_mort(self):
        #roll for death
        causes, risks = list(zip(*self.mort_risks))
        prob_nodeath = np.prod([1-p for p in risks])
        if random() > prob_nodeath:
            #roll for cause of death
            cause = causes[draw_dist(risks)]
            self.kill_patient(cause)
    def update_end_month(self):
        qol = self.inputs.qol
        
        #add costs for background care
        status = self.ss
        if status == SS_FORMER:
            months_quit = self.month - self.month_of_quit
            if months_quit < self.inputs.costs.bkgd_trans:
                status = SS_CURRENT
        self.add_cost("cost_bkgd",None,
                      self.inputs.costs.bkgd[self.gender,status,self.agecat], True)

        #add qol
        if self.update_qol:
            mult = qol.base[self.si, self.ss]
            if self.ss == SS_FORMER and (self.month - self.month_of_quit)<= qol.quit_duration[self.si]:
                mult *= qol.quit[self.si]
            self.add_qol(mult)

        #add outputs
        self.add_disc_out('lms',(self.ss,self.agecat, self.gender),1.0)
        self.add_disc_out('lms_SI',(self.si,),1.0)
        self.add_disc_out('qalms',(self.ss,self.agecat, self.gender),self.qol)


        for intv in range(NUM_INTERVENTIONS):
            if self.has_int[intv]:
                #add cost
                self.add_cost("cost_int",intv,self.inputs.costs.int_month[intv], True)
                #add output
                self.add_month_out("int_num_with",intv,1)

        for proph in range(NUM_PROPHS):
            if self.has_proph[proph]:
                #add cost
                self.add_cost("cost_proph",(proph,),self.inputs.costs.proph_month[proph], True)
                #add output
                self.add_month_out("proph_num_with",proph,1)
                
        for event in range(NUM_EVENTS):
            if self.events[event] == EVT_FULL:
                self.add_month_out('event_num_with',(event,self.ss,self.agecat,self.gender),1)
            else:
                self.add_month_out('event_num_without',event,1)

        #kill patient if max age reached
        if self.isalive and self.age / 12 >= self.inputs.sim.maxage:
            self.kill_patient(DTH_OLD_AGE)

        if self.isalive:
            self.add_month_out('num_alive',(self.ss,self.si), 1)
        else:
            self.add_month_out('num_deaths',(self.ss,self.si), 1)
        #advance age
        if self.isalive:
            self.month+=1
            self.age+=1

            #update disc factor
            self.disc_factor*=1/self.disc_mult_month

    def screen(self, testtype, event):
        #screen for pre events handles intial screen and schedules conf test
        pevt = self.inputs.pre_events

        #add cost
        self.add_cost("cost_screening",(event,),self.inputs.costs.screen[event], True)

        #add qol
        if self.update_qol:
            self.add_qol(self.inputs.qol.screen[event])
            
        true_status = self.events[event] == EVT_PRE
        if true_status:
            rate = pevt.screening_sensitivity
        else:
            rate = pevt.screening_specificity
            
        if random() <= rate[self.gender, event]:
            result = true_status
        else:
            result = not true_status

        self.add_month_out("pevent_screen_results", (event,int(true_status),int(result)),1)

        #prints to trace file if it exists
        if self.traces is not None:
            self.trace("\n**{0} Pre-event {1} Screening, Result: {2},Status: {2}".format(self.month,
                                                                                       self.event_names[event],
                                                                                       result, true_status))
        #schedule conf test
        if result:
            self.schedule_conf_test(event, self.month+pevt.screening_conf_delay[event])
            self.add_cost("cost_screening",event,self.inputs.costs.screen_pos[event], True)
        else:
            self.add_cost("cost_screening",event,self.inputs.costs.screen_neg[event], True)
    def conf_test(self,event):
        #handles a conf test
        pevt = self.inputs.pre_events
        self.month_of_conf[event] = None

        #add cost
        self.add_cost("cost_screening",event,self.inputs.costs.screen_conf[event], True)

        #add output
        self.add_month_out("pevent_conf",event,1)
        
        #add qol
        if self.update_qol:
            self.add_qol(self.inputs.qol.screen_conf[event])
        
        #add conf test mortality
        self.add_mort_risk(DTH_CONF_TEST,pevt.screening_conf_mort[event])


        #prints to trace file if it exists
        if self.traces is not None:
            self.trace("\n**{0} Pre-event {1} Confirmatory Test".format(self.month,self.event_names[event]))
        #change outcomes for true positives
        if self.events[event] == EVT_PRE:
            self.is_detected_pevent[event] = True
            #roll for revert to neutral
            if random() <= pevt.screening_outcome_neutral[event]:
                self.cure_pevent(event)
            #start proph if set
            if pevt.screening_outcome_proph[event]:
                proph_index = pevt.screening_outcome_proph[event] - 1
                if proph_index >= 0:
                    self.start_proph(proph_index)
            #start intervention if set
            if pevt.screening_outcome_intervention[event]:
                self.start_intervention(event)
            
    #add a mortality risk for that month
    def add_mort_risk(self, cause, prob, mult = 1):
        self.mort_risks.append((cause, prob*mult))
    #adds value to output
    def add_out(self, name, index, value, add_month = False):
        self.outputs.add_value(name, self.pid,index,value)

        #if true add to same name category in monthly costs
        if add_month:
            self.add_month_out(name, index, value)
            
    def add_list_out(self, name,value):
        self.outputs.add_list_value(name,value)
    def add_disc_out(self, name, index, value):
        out = self.disc_outputs
        out.add_value(name,self.pid,index, value, self.disc_factor)
    def add_month_out(self, name, index, value):
        self.curr_mo.add_value(name, index, value)
    def add_cost(self, name, index, value, add_month = False):
        out = self.disc_outputs
        disc_value = value*self.disc_factor
        out.add_value("overall_costs",self.pid,(self.ss, self.agecat, self.gender), value, self.disc_factor)
        out.add_value(name,self.pid,index, value, self.disc_factor)

        self.add_month_out("costs_disc", None, disc_value)

        #if true add to same name category in monthly costs
        if add_month:
            self.add_month_out(name, index, disc_value)
            
        
        
    def add_qol(self, qol_mult):
        self.qol*=qol_mult
    def set_proph_tox(self, proph):
        self.has_proph_tox_hist[proph] = True
        #add qol
        if self.update_qol:
            self.add_qol(self.inputs.qol.proph_tox[proph])
        #add output
        self.add_out("proph_num_tox",proph,1)

        #prints to trace file if it exists
        if self.traces is not None:
            self.trace("\n**{0} Toxicity Proph: {1} ".format(self.month,self.proph_names[proph]))
    def set_int_tox(self, intv):
        self.has_int_tox_hist[intv] = True
        #add qol
        if self.update_qol:
            self.add_qol(self.inputs.qol.int_tox[intv])

        #add output
        self.add_out("int_num_tox",intv,1)

        #prints to trace file if it exists
        if self.traces is not None:
            self.trace("\n**{0} Toxicity Intervention: {1} ".format(self.month,self.int_names[intv]))
    def schedule_conf_test(self, event, month_of_test):
        self.month_of_conf[event] = month_of_test

    def stop_proph(self, proph):
        self.has_proph[proph] = False
        #add output
        self.add_month_out("proph_stop",proph,1)

        #prints to trace file if it exists
        if self.traces is not None:
            self.trace("\n**{0} Stopping Proph: {1}".format(self.month, self.proph_names[proph]))
    def start_proph(self, proph, is_init = False):
        self.has_proph[proph] = True
        #add cost
        if not is_init:
            self.add_cost("cost_proph",(proph,),self.inputs.costs.proph_init[proph], True)
            #add output
            self.add_month_out("proph_start",proph,1)

            #prints to trace file if it exists
            if self.traces is not None:
                self.trace("\n**{0} Starting Proph: {1}".format(self.month, self.proph_names[proph]))
    def stop_intervention(self, intv):
        self.has_int[intv] = False

        #add output
        self.add_month_out("int_stop",intv,1)

        #prints to trace file if it exists
        if self.traces is not None:
            self.trace("\n**{0} Stopping Intervention: {1}".format(self.month, self.int_names[intv]))
    def start_intervention(self, intv, is_init = False):
        self.month_start_int[intv] = self.month
        self.has_int[intv] = True
        self.ever_start_any_int = True
        #add cost
        if not is_init:
            self.add_cost("cost_int",intv,self.inputs.costs.int_init[intv], True)
            #add output
            self.add_out("int_start",intv,1)
            self.add_month_out("int_start",intv,1)

            #prints to trace file if it exists
            if self.traces is not None:
                self.trace("\n**{0} Starting Intervention: {1}".format(self.month, self.int_names[intv]))
    def start_smoking(self):
        if self.ss ==SS_FORMER:
            self.add_list_out("quit_duration",self.month - self.month_of_quit)
        self.ss = SS_CURRENT
        
        #add output
        self.add_out("smoking_start",None,1, True)
        
        #prints to trace file if it exists
        if self.traces is not None:
            self.trace("\n**{0} Smoking Start".format(self.month))
    def quit_smoking(self):
        self.ss = SS_FORMER
        self.ever_quit= True
        self.month_of_quit = self.month
        self.agequit = self.age
        self.agequit_cat = get_age_cat(AGE_BRACKETS, self.agequit)

        #prints to trace file if it exists
        if self.traces is not None:
            self.trace("\n**{0} Smoking Quit".format(self.month))

        #check to see if should stop intervention
        for i in range(NUM_INTERVENTIONS):
            if self.has_int[i] and self.inputs.intervention.stop_on_quit[i]:
                self.stop_intervention(i)
                
        #add output
        self.add_out("smoking_quit",(self.agecat,self.gender),1, True)
        self.add_list_out("age_at_quit",self.agequit)
    def relapse_smoking(self):
        self.ss = SS_CURRENT
        self.ever_relapse = True
        #add output
        self.add_out("smoking_relapse",(self.agecat,self.gender),1, True)
        #prints to trace file if it exists
        if self.traces is not None:
            self.trace("\n**{0} Smoking Relapse".format(self.month))
    def get_event(self, event, is_prev=False):
        self.events[event] = EVT_FULL
        self.month_of_event[event] = self.month
        #add cost
        if not is_prev:
            self.add_cost("cost_event",event,self.inputs.costs.event_init[event], True)
            #add qol
            if self.update_qol:
                self.add_qol(self.inputs.qol.event_init[event])

            #prints to trace file if it exists
            if self.traces is not None:
                self.trace("\n**{0} Event {1}".format(self.month, self.event_names[event]))

        #add output
        if is_prev:
            self.add_out("event_num_prev",event,1)
        else:
            self.add_out("event_num_inc", event,1)
            self.add_month_out('event_inc',(event, self.ss, self.agecat, self.gender),1)
        
        #remove scheduled conf test
        self.month_of_conf[event] = None
        
    def cure_pevent(self, event):
        self.events[event] = EVT_NONE
        self.is_detected_pevent[event] = False
        #prints to trace file if it exists
        if self.traces is not None:
            self.trace("\n**{0} Pre-event {1} Cured".format(self.month,self.event_names[event]))
    def get_pevent(self, event):
        self.events[event] = EVT_PRE

        #add output
        self.add_out("pevent_num_inc", event,1)
        self.add_month_out('pevent_inc',(event, self.ss, self.agecat, self.gender),1)

        #prints to trace file if it exists
        if self.traces is not None:
            self.trace("\n**{0} Pre-event {1}".format(self.month,self.event_names[event]))
    def get_event_comp(self, event):
        self.month_of_comp[event] = self.month
        #add cost
        self.add_cost("cost_comp",(event,),self.inputs.costs.comp[event], True)

        #add qol
        if self.update_qol:
            self.add_qol(self.inputs.qol.comp[event])

        #add output
        self.add_out("event_comp", event,1)
        self.add_month_out('event_comp',(event, self.ss, self.agecat, self.gender),1)

        #prints to trace file if it exists
        if self.traces is not None:
            self.trace("\n**{0} Event {1} Complication".format(self.month,self.event_names[event]))
    #run single time step(month)
    def run_step(self):
        self.update_start_month()
        self.update_interventions()
        self.update_smoke_start()
        self.update_smoke_quit()
        self.update_smoke_relapse()
        self.update_prophs()
        self.update_pre_events()
        self.update_pevent_screening()
        self.update_events()
        self.update_nathist()
        self.update_mort()
        self.update_end_month()
        
    def run_until_death(self):
        while(self.isalive):
            self.run_step()

        if self.traces is not None:
            self.traces[self.pid] = self.trace_text

    def kill_patient(self,cause):
        self.isalive = False
        self.causeofdeath = cause

        #record output
        self.add_out("death_dist",(self.ss,self.si),1)
        for event in range(NUM_EVENTS):
            if self.events[event] == EVT_FULL:
                self.add_out("event_mths",event,self.month - self.month_of_event[event])
        self.add_out("death_causes",cause,1, True)

        if self.ever_quit:
            self.add_out("smoking_ever_quit",None, 1)
        if self.ever_relapse:
            self.add_out("smoking_ever_relapse",None, 1)
        
        if self.ever_start_any_int:  
            self.add_out("int_start",None,1)
            
        if self.traces is not None:
            self.trace("\n**{0} DEATH {1}".format(self.month,self.death_names[cause]))
            self.trace("\nEND PATIENT {0}".format(self.pid), False)
    def trace(self, text, include_lm = True):
        #prints to trace file
        out = self.disc_outputs.disc
        if include_lm:
            text+=", LM {0:.2f}, QA {1:.2f}, $ {2:.2f}".format(out['lms'][self.pid].sum(),
                                                   out['qalms'][self.pid].sum(),
                                                   out['overall_costs'][self.pid].sum())

        self.trace_text+=text
        #self.traces.write(text)
        
    def __repr__(self):
        s = "Patient {0}".format(self.pid)
        s += "\nGender:{0}".format(GENDER_STRS[self.gender])
        s += "\nAge:{0}({1}Y,{2}M)".format(self.age,*convertage(self.age))
        s += "\nSmoking:{0}, {1}".format(SS_STRS[self.ss],
                                                SI_STRS[self.si])
        return s
if __name__ == "__main__":
    excelfile = "Smoking Model Inputs.xlsm"
    runsize = 100
    i = Inputs()
    i.load_excel(excelfile)
    o = Outputs(runsize, i)
    
    do = DiscOutputs(runsize, i)
    mo = MonthlyOutputs()
    num = 0
    tracepath = "test.smtrace"
    with open(tracepath, 'w') as ftrace:
        for n in range(runsize):
            if n%100==0:
                print(n)
            
        
            if n <50:
                traces = ftrace
            else:
                traces = None
            
            p = Patient(n, i, o, do, mo, traces)
            p.run_until_death()
        
    write_output("test.smout",o, do, mo)        
