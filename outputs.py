"""
Handles reading and manipulating outputs for the smoking model

@author: Taige Hou(thou1@partners.org)
"""

from enums import *
import numpy as np

######################################################################
#(varname, shape, init,type)
OUTPUTS = [("init_dist_smoking",[len(SS),len(SI)],0,int),
           ("init_dist_gender",[len(GENDERS)],0,int),
           ("init_age",[],0,int),
           ("death_dist",[len(SS),len(SI)],0,int),
           ("death_causes",[len(DTH_CAUSES)],0,int),
           ("event_num_prev",[NUM_EVENTS],0,int),
           ("event_num_inc",[NUM_EVENTS],0,int),
           ("event_mths",[NUM_EVENTS], 0, int),
           ("pevent_num_inc",[NUM_EVENTS],0,int),
           ("pevent_mths",[NUM_EVENTS], 0, int),
           ("event_comp", [NUM_EVENTS],0,int),

            ("int_num_tox", [NUM_INTERVENTIONS], 0, int),
           ("int_start",[NUM_INTERVENTIONS], 0, int),
           
           ("proph_num_tox", [NUM_PROPHS], 0, int),

           ("smoking_start",[],0, int),
           ("smoking_quit",[len(AGE_BRACKETS),len(GENDERS)],0,int),
           ("smoking_relapse", [len(AGE_BRACKETS),len(GENDERS)],0,int),
           ("smoking_ever_quit",[],0,int),
           ("smoking_ever_relapse",[],0,int),
           ("int_ever_start",[], 0, int),
            ("quit_duration",[],0,list),
           ("age_at_quit",[],0,list),
           ]

DISC_OUTPUTS =[("lms", [len(SS),len(AGE_BRACKETS), len(GENDERS)],0.0, float),
               ("lms_SI",[len(SI)],0.0,float),
                ("qalms", [len(SS),len(AGE_BRACKETS), len(GENDERS)],0.0, float),
                ("cost_event", [NUM_EVENTS,],0.0, float),
               ("cost_comp", [NUM_EVENTS,], 0.0, float),
               ("cost_proph", [NUM_PROPHS,],0.0, float),
               ("cost_int", [NUM_INTERVENTIONS,], 0.0, float),
               ("cost_bkgd", [],0.0,float),
               ("cost_screening",[NUM_EVENTS],0.0,float),

           ]

MONTH_OUTPUTS = [("num_alive",(len(SS),len(SI)),0,int),
                 ("num_deaths",(len(SS),len(SI)),0,int),
                 ("death_causes",[len(DTH_CAUSES)],0,int),
                 ("costs_disc",[],0.0,float),

                 ("event_num_with",[NUM_EVENTS, len(SS),len(AGE_BRACKETS),len(GENDERS)],0,int),
                 ("event_num_without",[NUM_EVENTS],0,int),
                 ("event_inc",[NUM_EVENTS, len(SS),len(AGE_BRACKETS),len(GENDERS)],0,int),
                 ("pevent_inc",[NUM_EVENTS, len(SS),len(AGE_BRACKETS),len(GENDERS)],0,int),
                 ("event_comp",[NUM_EVENTS, len(SS),len(AGE_BRACKETS),len(GENDERS)],0,int),
                 ("pevent_screen_results",[NUM_EVENTS,2,2],0,int), #(truth, obsv)
                 ("pevent_conf",[NUM_EVENTS],0,int),

                 ("smoking_start",[],0, int),
                 ("smoking_quit",[len(AGE_BRACKETS),len(GENDERS)],0,int),
                 ("smoking_relapse", [len(AGE_BRACKETS),len(GENDERS)],0,int),
                 
                 ("cost_event", [NUM_EVENTS,],0.0, float),
                 ("cost_comp", [NUM_EVENTS,], 0.0, float),
                 ("cost_screening",[NUM_EVENTS],0.0,float),
                 ("cost_proph", [NUM_PROPHS,],0.0, float),
                 ("cost_int", [NUM_INTERVENTIONS,], 0.0, float),
                 ("cost_bkgd", [],0.0,float),

                 ("proph_num_with",[NUM_PROPHS],0,int),
                 ("proph_start",[NUM_PROPHS],0,int),
                 ("proph_stop",[NUM_PROPHS],0,int),

                 ("int_num_with",[NUM_INTERVENTIONS],0,int),
                 ("int_start",[NUM_INTERVENTIONS],0,int),
                 ("int_stop",[NUM_INTERVENTIONS],0,int),

                 
    ]
    
#Overall Outputs (undiscounted)
class Outputs(dict):
    def __init__(self, num_patients, inputs, outputs = OUTPUTS):
        super(Outputs,self).__init__()
        self.num_patients = num_patients
        self.init_outputs(outputs)
        self.inputs = inputs
        self.disc_rate = inputs.sim.disc_rate_year
    def init_outputs(self, outputs):
        for output in outputs:
            self.add_output(output)
    def add_output(self, params):
        varname, shape, value, typ =params
        if typ == list:
            #dont create numpy array this will be a dynamically grown list
            self[varname] = []
        else:
            new_shape = [self.num_patients]+shape
            self[varname] = np.full(new_shape, value, typ)
    def add_value(self, *args):
        varname, pid, index, value = args
        if index is not None:
            self[varname][pid][index]+=value
        else:
            self[varname][pid]+=value
    def add_list_value(self, *args):
        name, value = args
        self[name].append(value)
    def __getattr__(self, name):
        return self[name]    
    def __setattr__(self,name, value):
        self[name] = value


    
#stores both discounted and undiscounted outputs
class DiscOutputs(object):
    def __init__(self, num_patients, inputs):
        self.disc = Outputs(num_patients, inputs, DISC_OUTPUTS)
        self.undisc = Outputs(num_patients, inputs, DISC_OUTPUTS)
        self.disc.add_output(("overall_costs",[len(SS),len(AGE_BRACKETS), len(GENDERS)],0.0, float))
        self.undisc.add_output(("overall_costs",[len(SS),len(AGE_BRACKETS), len(GENDERS)],0.0, float))
    def add_value(self, *args):
        varname, pid, index, value, disc_factor = args
        if index is not None:
            self.undisc[varname][pid][index]+=value
            self.disc[varname][pid][index]+=value*disc_factor
        else:
            self.undisc[varname][pid]+=value
            self.disc[varname][pid]+=value*disc_factor



#class so store outputs for a single month
class MonthOut(dict):
    def __init__(self, outputs = MONTH_OUTPUTS):
        super(MonthOut,self).__init__()
        for output in outputs:
            self.add_output(output)
    def add_output(self, params):
        varname, shape, value, typ =params
        self[varname] = np.full(shape, value, typ)
    def add_value(self, varname, index, value):
        if index is not None:
            self[varname][index]+=value
        else:
            self[varname]+=value
    def __getattr__(self, name):
        return self[name]    
    def __setattr__(self,name, value):
        self[name] = value

#class to store all monthly outputs
class MonthlyOutputs(list):
    def __init__(self):
        super(MonthlyOutputs,self).__init__()
    def add_month(self):
        mo = MonthOut()
        self.append(mo)
        return mo
    
############################################################################################################
def write_avg(fout, label, a, axis):
    a = a.sum(axis =axis)
    mean = a.mean()
    std = a.std()

    fout.write("\n\t{0}\t{1}\t{2}".format(label, mean,std))

def write_single(fout, header, a, avg=True, ismonth = False):
    if not ismonth:
        if avg:
            a = a.mean(axis=0)
        else:
            a= a.sum(axis=0)

    fout.write("\n\t{0}\t{1}".format(header, a))
               
def write_1d_array(fout, header, labels, a, avg= True, write_labels = True, ismonth=False):
    if not ismonth:
        #avg across patients
        if avg:
            a = a.mean(axis=0)
        else:
            a=a.sum(axis=0)

    fout.write("\n")
    if write_labels:
        fout.write("\n\t\t"+"\t".join(labels))
    fout.write("\n\t{0}".format(header))
    for x in range(a.size):
        fout.write("\t{0}".format(a[x]))
            
def write_2d_array(fout, header, labels, a, order, avg= True, ismonth=False):
    if not ismonth:
        #avg across patients
        if avg:
            a = a.mean(axis=0)
        else:
            a=a.sum(axis=0)
    a= a.transpose(order)
    i,j=order
    dims = a.shape
    fout.write("\n")
    fout.write("\n\t{0}\t".format(header)+"\t".join(labels[j]))
    for x in range(dims[0]):
        fout.write("\n\t{0}".format(labels[i][x]))
        for y in range(dims[1]):
            fout.write("\t{0}".format(a[x,y]))
def write_3d_array(fout, header, labels, a, order, avg= True, ismonth = False):
    if not ismonth:
        #avg across patients
        if avg:
            a = a.mean(axis=0)
        else:
            a=a.sum(axis=0)
    a= a.transpose(order)
    i,j,k=order
    dims = a.shape
    fout.write("\n")
    fout.write("\n\t{0}\t\t".format(header)+"\t".join(labels[k]))
    for x in range(dims[0]):
        fout.write("\n\t{0}".format(labels[i][x]))
        for y in range(dims[1]):
            if y == 0:
                fout.write("\t{0}".format(labels[j][y]))
            else:
                fout.write("\n\t\t{0}".format(labels[j][y]))
            for z in range(dims[2]):
                fout.write("\t{0}".format(a[x,y,z]))
                
   
#writes output file
def write_output(filepath,out,dout, mout):
    #write header
    with open(filepath,'w') as fout:
        fout.write("POPULATION SUMMARY MEASURES")
        fout.write("\n\tRun Size\t{0}".format(out.num_patients))
        fout.write("\n\tDisc Rate\t{0}".format(dout.disc.disc_rate))
        
        fout.write("\n\n\tOutcome\tMean\tStd Dev")
        write_avg(fout, "Costs Disc",dout.disc.overall_costs,(1,2,3))
        write_avg(fout, "Life Months Disc",dout.disc.lms,(1,2,3))
        write_avg(fout, "Quality-Adjusted Life Months Disc",dout.disc.qalms,(1,2,3))
        fout.write("\n")
        write_avg(fout, "Costs Undisc",dout.undisc.overall_costs,(1,2,3))
        write_avg(fout, "Life Months Undisc",dout.undisc.lms,(1,2,3))
        write_avg(fout, "Quality-Adjusted Life Months Undisc",dout.undisc.qalms,(1,2,3))

        write_3d_array(fout, "Costs",(SS_STRS,AGE_BRACKET_STRS,GENDER_STRS),dout.disc.overall_costs,(2,0,1))
        write_3d_array(fout, "Life Months",(SS_STRS,AGE_BRACKET_STRS,GENDER_STRS),dout.disc.lms,(2,0,1))
        write_3d_array(fout, "Quality-Adjusted Life Months",(SS_STRS,AGE_BRACKET_STRS,GENDER_STRS),dout.disc.qalms,(2,0,1))

        write_1d_array(fout, "Life Months",(SI_STRS), dout.disc.lms_SI)
        
        fout.write("\nINITIAL DISTRIBUTIONS")
        write_2d_array(fout, "Smoking",(SS_STRS,SI_STRS),out.init_dist_smoking,(0,1), False)
        write_1d_array(fout, "Gender",(GENDER_STRS),out.init_dist_gender, False)

        fout.write("\n\t\tMean\tStd Dev")
        write_avg(fout,"Age",out.init_age,())

        
        write_deaths(fout, out, dout)
        write_events(fout,out,dout)
        write_smoking(fout,out,dout)
        write_interventions(fout, out, dout)
        write_prophs(fout, out, dout)
        write_costs(fout, out, dout)
        write_months(fout, out, mout)
def write_deaths(fout, out, dout):
    fout.write("\nDEATH DISTRIBUTIONS")
    write_2d_array(fout, "Death",(SS_STRS, SI_STRS),out.death_dist, (0,1), False)
    death_causes = ["Nat Hist","Old Age", "Proph Tox", "Intervention Tox", "Confirmatory Test",]
    death_causes.extend(["{0}".format(out.inputs.sim.event_names[i]) for i in range(NUM_EVENTS)])
    death_causes.extend(["{0} Complication".format(out.inputs.sim.event_names[i]) for i in range(NUM_EVENTS)])

    write_1d_array(fout, "Causes of Death", death_causes, out.death_causes,False)
def write_events(fout, out, dout):
    fout.write("\nEVENTS")

    #headers
    headers = ("mean months with event","num prevalent","num incident",
               "mean months with pre-event","num incident pre-event",
               "num complications",)
    fout.write("\n\t\t"+"\t".join(headers))
    values = (np.ma.masked_equal(out.event_mths,0).mean(axis=0),
              out.event_num_prev.sum(axis=0),
              out.event_num_inc.sum(axis=0),
              np.ma.masked_equal(out.pevent_mths,0).mean(axis=0),
              out.pevent_num_inc.sum(axis=0),
              out.event_comp.sum(axis=0),
              )

    for event in range(NUM_EVENTS):
        fout.write("\n\t{0}".format(out.inputs.sim.event_names[event]))
        for value in values:
            fout.write("\t{0}".format(value[event]))

def write_smoking(fout, out, dout):
    fout.write("\nSMOKING")
    
    write_single(fout, "Num Start",out.smoking_start, False)
    write_2d_array(fout, "Quit Events",(AGE_BRACKET_STRS, GENDER_STRS),out.smoking_quit, (1,0), False)
    write_2d_array(fout, "Relapse Events", (AGE_BRACKET_STRS, GENDER_STRS), out.smoking_relapse, (1, 0), False)
    write_single(fout, "Ever Quit",out.smoking_ever_quit, False)
    write_single(fout, "Ever Relapse",out.smoking_ever_relapse, False)
    write_single(fout, "Ever Intervention",out.int_ever_start, False)

    if out.quit_duration:
        q3,q1 = np.percentile(out.quit_duration,[75,25])
        iqr = q3-q1
        values = (np.mean(out.quit_duration),
                  np.std(out.quit_duration),
                  min(out.quit_duration),
                  max(out.quit_duration),
                  iqr)
    else:
        values = ['--']*5
    fout.write("\n\t\tMean\tStd Dev\tMin\tMax\tIQR")
    fout.write("\n\tQuit Duration (Of those who relpase)\t{0}\t{1}\t{2}\t{3}\t{4}".format(*values))
    if out.age_at_quit:
        values = (np.mean(out.age_at_quit),
                    np.std(out.age_at_quit))
    else:
        values = ['--']*2
    fout.write("\n\tAge at Quit\t{0}\t{1}".format(*values))
    
    
def write_interventions(fout, out, dout):
    fout.write("\nINTERVENTIONS")

    headers = ("num start events", "tox events", "intervention costs")
    values = (out.int_start.sum(axis=0),
            out.int_num_tox.sum(axis=0),
              dout.disc.cost_int.sum(axis=0),
        )
    fout.write("\n\t\t"+"\t".join(headers))
    for intv in range(NUM_INTERVENTIONS):
        fout.write("\n\t{0}".format(out.inputs.intervention.int_names[intv]))
        for value in values:
            fout.write("\t{0}".format(value[intv]))

def write_prophs(fout, out, dout):
    fout.write("\nPROPHS")

    headers = ("tox events", "proph costs")
    values = (out.proph_num_tox.sum(axis=0),
              dout.disc.cost_proph.sum(axis=0),
        )
    fout.write("\n\t\t"+"\t".join(headers))
    for proph in range(NUM_PROPHS):
        fout.write("\n\t{0}".format(out.inputs.prophs.proph_names[proph]))
        for value in values:
            fout.write("\t{0}".format(value[proph]))

def write_costs(fout, out, dout):
    fout.write("\nCOSTS")

    write_1d_array(fout, "Event Costs",out.inputs.sim.event_names, dout.disc.cost_event, False, True)
    write_1d_array(fout, "Event Complication Costs",out.inputs.sim.event_names, dout.disc.cost_comp, False, False)
    write_1d_array(fout, "Event Screening Costs",out.inputs.sim.event_names, dout.disc.cost_screening, False, False)

    write_single(fout, "Background Care Costs",dout.disc.cost_bkgd, False)

def write_months(fout, out, mout):
    for n,mo in enumerate(mout):
        fout.write("\nCohort Summary for Month {0}".format(n))
        write_single(fout, "Num Alive", mo.num_alive.sum(), ismonth = True)
        write_2d_array(fout, "Num Alive",(SS_STRS,SI_STRS), mo.num_alive, (0,1), ismonth=True)
        fout.write("\n")
        write_single(fout, "Num Deaths", mo.num_deaths.sum(), ismonth = True)
        write_2d_array(fout, "Num Deaths",(SS_STRS,SI_STRS), mo.num_deaths, (0,1), ismonth=True)

        death_causes = ["Nat Hist","Old Age", "Proph Tox", "Intervention Tox", "Confirmatory Test",]
        death_causes.extend(["{0}".format(out.inputs.sim.event_names[i]) for i in range(NUM_EVENTS)])
        death_causes.extend(["{0} Complication".format(out.inputs.sim.event_names[i]) for i in range(NUM_EVENTS)])
        write_1d_array(fout, "Causes of Death", death_causes, mo.death_causes,ismonth = True)
    
        fout.write("\n")
        write_single(fout, "Monthly Costs", mo.costs_disc, ismonth = True)
        write_single(fout, "Background Care Costs", mo.cost_bkgd, ismonth = True)
        write_single(fout, "Num Start", mo.smoking_start, ismonth = True)
        write_2d_array(fout, "Quit Events", (AGE_BRACKET_STRS, GENDER_STRS), mo.smoking_quit, (1, 0), ismonth=True)
        write_2d_array(fout, "Relapse Events", (AGE_BRACKET_STRS, GENDER_STRS), mo.smoking_relapse, (1, 0), ismonth=True)

        #events
        for event in range(NUM_EVENTS):
            fout.write("\n")
            write_single(fout,"Num Alive without {0}".format(out.inputs.sim.event_names[event]),
                          mo.event_num_without[event],ismonth = True)
            write_3d_array(fout, "Num Alive with {0}".format(out.inputs.sim.event_names[event]),
                           (SS_STRS,AGE_BRACKET_STRS,GENDER_STRS), mo.event_num_with[event], (2,0,1), ismonth = True)
            write_3d_array(fout, "Num Incident {0}".format(out.inputs.sim.event_names[event]),
                           (SS_STRS,AGE_BRACKET_STRS,GENDER_STRS), mo.event_inc[event], (2,0,1), ismonth = True)
            write_3d_array(fout, "Num Incident Pre-event {0}".format(out.inputs.sim.event_names[event]),
                           (SS_STRS,AGE_BRACKET_STRS,GENDER_STRS), mo.pevent_inc[event], (2,0,1), ismonth = True)

            write_2d_array(fout, "Pre-event Screening Results {0}".format(out.inputs.sim.event_names[event]),
                           (("True Neg","True Pos"),("Obsv Neg","Obsv Pos")),mo.pevent_screen_results[event],(0,1),ismonth=True)

            write_single(fout,"Pre-event Confirmatory Tests {0}".format(out.inputs.sim.event_names[event]),
                         mo.pevent_conf[event],ismonth = True)

            write_3d_array(fout, "Num Complications {0}".format(out.inputs.sim.event_names[event]),
                           (SS_STRS,AGE_BRACKET_STRS,GENDER_STRS), mo.event_comp[event], (2,0,1), ismonth = True)
            write_single(fout,"Cost Events {0}".format(out.inputs.sim.event_names[event]),
                         mo.cost_event[event],ismonth = True)
            write_single(fout,"Cost Complications {0}".format(out.inputs.sim.event_names[event]),
                         mo.cost_comp[event],ismonth = True)
            write_single(fout,"Cost Pre-event Screening {0}".format(out.inputs.sim.event_names[event]),
                         mo.cost_screening[event],ismonth = True)
            
        #Interventions
        fout.write("\n")
        headers = ("num on intervention", "num start intervention","num stop intervention", "intervention cost")
        values = (mo.int_num_with,
                  mo.int_start,
                  mo.int_stop,
                  mo.cost_int,
                )
        fout.write("\n\t\t"+"\t".join(headers))
        for intv in range(NUM_INTERVENTIONS):
            fout.write("\n\t{0}".format(out.inputs.intervention.int_names[intv]))
            for value in values:
                fout.write("\t{0}".format(value[intv]))
                

        #Prophs
        fout.write("\n")
        headers = ("num on proph", "num start proph","num stop proph", "proph cost")
        values = (mo.proph_num_with,
                  mo.proph_start,
                  mo.proph_stop,
                  mo.cost_proph,
                )
        fout.write("\n\t\t"+"\t".join(headers))
        for proph in range(NUM_PROPHS):
            fout.write("\n\t{0}".format(out.inputs.prophs.proph_names[proph]))
            for value in values:
                fout.write("\t{0}".format(value[proph]))
if __name__=="__main__":
    from inputs import *
    excelfile = "../Smoking Model Inputs.xlsm"
    i = Inputs()
    i.load_excel(excelfile)
    o = Outputs(1000, i)
    do = DiscOutputs(1000, i)
    mo = MonthlyOutputs()
    mo.add_month()
    costs = do.disc.overall_costs
    write_output("test.smout",o, do, mo)

    
    
