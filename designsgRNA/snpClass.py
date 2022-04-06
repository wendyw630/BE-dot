### snp in snp_class containing ReadingFrame
class snp_define:
    def __init__(self,jobid,seq5,seq3,wt,mut,readingframe):
        self.jobid=jobid
        self.seq5=seq5[25:]
        self.seq3=seq3[:25]
        self.mut=mut
        self.wt=wt
        self.readingframe=readingframe
    def __repr__(self):
        return "Id's:" + self.jobid + " , " +\
               " \nsequence:" + self.seq5 + self.mut + self.seq3 +"\nwt:" + self.wt

class supo_snp_define:
    ## wt,mut -> wt,edited
    def __init__(self,jobid,seq5,seq3,wt,readingframe):
        self.jobid=jobid
        self.seq5=seq5
        self.seq3=seq3
        self.wt=wt
        self.readingframe=readingframe
    def __repr__(self):
        return "Id's:" + self.jobid + " , " +\
               " \nsequence:" + self.seq5 +" "+ self.seq3 +"\nwt:" + self.wt