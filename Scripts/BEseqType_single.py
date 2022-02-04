#AUTHOR: Shiri Almog , shirialmog1@gmail.com
class BEseqType:
    def __init__(self,jobID,seq5,seq3, wt,mutation,readingFrame): # wzl_change: del(snpID, aaPosition,geneName,geneID), add(jobID)
        self.jobID = jobID
        self.seq5=seq5
        self.seq3=seq3
        self.wt = wt
        self.mutation=mutation
        self.readingFrame = readingFrame


    def __repr__(self):
        return "Id's:"+self.jobID+" , "\
               " \nsequence:"+self.seq5+ self.mutation + self.seq3+\
               " \nwt:"+self.wt
