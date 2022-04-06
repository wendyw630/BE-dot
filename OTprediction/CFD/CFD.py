# === SOURCE CODE cfd-score-calculator.py provided by John Doench =====
import pickle, re
from os.path import basename, join, splitext, isfile, dirname


def get_mm_pam_scores():

    # dataDir = join(dirname(__file__), 'CFD_Scoring')
    dataDir = r'D:\00_project\predict\BE_dot_02\OT_predict\CFD' ## wzl_note: the orientation of /
    mm_scores = pickle.load(open(join(dataDir, 'mismatch_score.pkl'), 'rb'))
    pam_scores = pickle.load(open(join(dataDir, 'pam_score.pkl'), 'rb'))
    return (mm_scores, pam_scores)
    # except:
    # raise Exception("Could not find file with mismatch scores or PAM scores")

# Reverse complements a given string
def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)

# Calculates CFD score
def calc_cfd(wt, sg): ### wzl_changed, orig: (wt,sg,pam)
    '''
    :param wt:
    :param sg:
    :return:
    '''
    """
    >>> calc_cfd("GGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGGGG", "GG")
    1.0
    """
    mm_scores,pam_scores = get_mm_pam_scores()
    score = 1
    sg = sg.replace('T', 'U')
    wt = wt.replace('T', 'U') ## wzl_note
    s_list = list(sg)
    wt_list = list(wt)
    # print mm_scores
    for i, sl in enumerate(s_list):
        # print i, sl, wt_list[i], mm_scores
        if wt_list[i] == sl:
            score *= 1
        else:
            key = 'r' + wt_list[i] + ':d' + revcom(sl) + ',' + str(i + 1)
            score *= mm_scores[key]
    #score *= pam_scores[pam] ### wzl_del
    return (score)

mm_scores, pam_scores = None, None
def calcCfdScore(guideSeq, otSeq):

    """ 
    based on source code provided by John Doench
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGAAAGGG")
    0.4635989007074176
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGGGGGGG")
    1.0
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "aaaaGaGaGGGGGGGGGGGGGGG")
    0.5140384614450001
    """
    global mm_scores  ### wzl_changed: orignal: global mm_scores, pam_scores
    if mm_scores is None:
        mm_scores = get_mm_pam_scores() ### wzl_changed, orig: mm_scores, pam_scores =
    wt = guideSeq.upper()[:20]
    off = otSeq.upper()[:20]
    m_wt = re.search('[^ATCG]', wt) ## wzl_note: [^AGCT]: not AorGorCorT
    m_off = re.search('[^ATCG]', off)
    if (m_wt is None) and (m_off is None):
        #pam = off[-2:] ## wzl_del #set len(pam)=3 by default, ............
        #sg = off[:-3] ## wzl_note: wt=guideSeq==23nt;   off=otSeq==23nt,off=sg+pam
        cfd_score = calc_cfd(wt, off) ### wzl_changed
        return cfd_score
        # print "CFD score: "+str(cfd_score)


"""
guideSeq='AAATCTTAGAGCGTCCCATC'
otPamSeq='AAATCTTAGAGcGTCCCATCTGG'
cfd_score = calcCfdScore(guideSeq, otPamSeq)
print("cfd_score",cfd_score)
"""


