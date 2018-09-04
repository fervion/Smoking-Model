"""
The main simulation object for smoking model

@author: Taige Hou(thou1@partners.org)
"""

from inputs import *
from outputs import *
from enums import *
from patient import *
import os, sys,threading
from glob import glob
import random

######################################################################
#patient thread
class PatientThread(threading.Thread):
    def __init__(self,tID,pargs):
        threading.Thread.__init__(self)
        self.tID = tID
        self.pargs = pargs
    def run(self):
        p = Patient(*self.pargs)
        p.run_until_death()

######################################################################
#Main simulation object
class Sim(object):
    def __init__(self):
        self.inputs = None
        self.traces = {}
    def load_inputs_xl(self, filepath):
        self.input_path = filepath

        self.inputs = Inputs()
        self.inputs.load_excel(filepath)

        self.runsize = self.inputs.sim.runsize
        self.is_rand_seed = self.inputs.sim.rand_seed

        self.outputs = Outputs(self.runsize, self.inputs)

        self.disc_outputs = DiscOutputs(self.runsize, self.inputs)

        self.monthly_outputs = MonthlyOutputs()

    def load_inputs_text(self, filepath):
        self.input_path = filepath

        self.inputs = Inputs()
        self.inputs.load_txt(filepath)

        self.runsize = self.inputs.sim.runsize
        self.is_rand_seed = self.inputs.sim.rand_seed

        self.outputs = Outputs(self.runsize, self.inputs)

        self.disc_outputs = DiscOutputs(self.runsize, self.inputs)

        self.monthly_outputs = MonthlyOutputs()
    def convert_excel_to_text(self,excelfile,textfile):
        self.load_inputs_xl(excelfile)
        self.inputs.save_txt(textfile)
    #main loop
    def run(self):
        threads = []
        tracepath = os.path.splitext(self.input_path)[0]+".smtrace"
        for n in range(self.runsize):

            if not self.is_rand_seed:
                random.seed(n)
            if n%100==0:
                print(n)
            if n < 50:
                traces = self.traces
            else:
                traces = None

            t = PatientThread(n,[n,self.inputs,self.outputs, self.disc_outputs, self.monthly_outputs, traces])
            t.start()
            threads.append(t)
            #p = Patient(n,self.inputs,self.outputs, self.disc_outputs, self.monthly_outputs, trace_file)
            #p.run_until_death()

        for t in threads:
            t.join()

        with open(tracepath, 'w') as ftrace:
            for pid in sorted(self.traces):
                ftrace.write(self.traces[pid])

        #write output
        outpath = os.path.splitext(self.input_path)[0]+".smout"
        write_output(outpath,self.outputs, self.disc_outputs, self.monthly_outputs)
if __name__ == "__main__":
    s = Sim()
    if len(sys.argv)>= 2 and sys.argv[1]=="text":
        #convert excel file to text
        print(sys.argv)
        excelfile = []
        textfile = []
        i = 2
        while True:
            excelfile.append(sys.argv[i])
            if sys.argv[i].endswith(".xlsm"):
                break
            i+=1
        i+=1
        while True:
            textfile.append(sys.argv[i])
            if i==len(sys.argv)-1:
                break
            i+=1
        excelfile = " ".join(excelfile)
        textfile= " ".join(textfile)
        s.convert_excel_to_text(excelfile,textfile)
    else:
        if len(sys.argv) ==1:
            excelfiles = glob("*.xlsm")
            textfiles = glob("*.smin")
        else:
            folder = sys.argv[1]
            excelfiles = glob(os.path.join(folder,"*.xlsm"))
            textfiles = glob(os.path.join(folder,"*.smin"))

        for excelfile in excelfiles:
            print(excelfile)
            try:
                s.load_inputs_xl(excelfile)
            except:
                pass
            else:
                s.run()

        for textfile in textfiles:
            print(textfile)
            try:
                s.load_inputs_text(textfile)
            except IndexError:
                pass
            else:
                s.run()