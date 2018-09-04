from random import random
from scipy.stats import truncnorm
import numpy as np
import re,operator

NUM_EVENTS = 5
NUM_PROPHS = 10
NUM_INTERVENTIONS = 10

GENDERS = FEMALE, MALE = (0,1)
GENDER_STRS = ("Female","Male")

SS = SS_NEVER, SS_FORMER, SS_CURRENT = (0,1,2)
SS_STRS = ("Never","Former","Curent")
SI = SI_LIGHT, SI_MODERATE, SI_HEAVY = (0,1,2)
SI_STRS = ("Light","Moderate","Heavy")

SCREEN_REG, SCREEN_BACK = (0,1)

EVENT_STAGES = EVT_NONE,EVT_PRE,EVT_FULL = (0,1,2)
DTH_NAMES = ["DTH_NAT_HIST", "DTH_OLD_AGE", "DTH_TOX_PROPH", "DTH_TOX_INT","DTH_CONF_TEST"]
DTH_NAMES.extend(["DTH_EVENT_{0}".format(i) for i in range(NUM_EVENTS)])
DTH_NAMES.extend(["DTH_EVENT_COMP_{0}".format(i) for i in range(NUM_EVENTS)])

for i,varname in enumerate(DTH_NAMES):
    globals()[varname] = i
DTH_CAUSES = tuple(range(len(DTH_NAMES)))


AGE_BRACKETS = (20,25,30,35,40,
                45,50,55,60,65,
                70,75,80,85,np.inf)
AGE_BRACKET_STRS = ("<20","20-24","25-29","30-34",
                    "35-39","40-44","45-49","50-54",
                    "55-59","60-64","65-69","70-74",
                    "75-79","80-84","85+")



#draws from a weighted distribution where distribution is a dictionary
#with choices to weights.  or a list of probs. weights are not normalized
def draw_dist(dist):
    if isinstance(dist, dict):
        sumWeights = sum(dist.values())
    else:
        sumWeights = sum(dist)        
        
    randNum = random()*sumWeights
    culNum = 0

    if isinstance(dist, dict):
        items = iter(dist.items())
    else:
        items = enumerate(dist)
        
    for key, value in items:
        culNum += value
        if randNum <= culNum:
            return key

#draws from truncated normal dist
def draw_trunc_norm(mean, std, low, high):
    if std == 0:
        return mean
    if low == -1:
        low == -np.inf
    if high == -1:
        high = np.inf
    a, b = (low - mean)/float(std), (high - mean)/float(std)
    return truncnorm(a,b, loc=mean, scale=std).rvs()

def convertage(agemonths):
    return divmod(agemonths,12)

#gets the index of the age category from the bracket
#bracket format can be years y or months m
def get_age_cat(bracket, agemonths, ageformat = "y"):
    if ageformat == "y":
        age = agemonths/12
    else:
        age = agemonths
    for i,upp in enumerate(bracket):
        if age < upp:
            return i
    return i

#gets age category for tripartition e.g. age<=N1, N1<age<=N2, age>N2
def get_tribound_cat(bounds, agemonths, ageformat = "y", boundsformat = ("<=","<=")):
    if ageformat == "y":
        age = agemonths/12
    else:
        age = agemonths

    ops = {'>': operator.gt,
           '<': operator.lt,
           '>=': operator.ge,
           '<=': operator.le,
           '=': operator.eq,
           }
    for i,op in enumerate(boundsformat):
        if ops[op](age, bounds[i]):
            return i

    return i+1
if __name__ == "__main__":
    import numpy as np
    print(draw_trunc_norm(20,1,-1,-1))
    
