#/usr/bin/python3


import Bio.SeqIO
import re
import sys

if (len(sys.argv)) < 4:
    print("Usage:", sys.argv[0], "cosmo.stdout1 cosmo.stdout2 flank_cmp_size")
    sys.exit(1)

    
flankid_re = re.compile(r"var_(\d+)_(\dp)_flank")
branchid_re = re.compile(r"var_(\d+)_branch_(\d)")
def parse_cortex(cortex_fname):
    """Parse a cortex bubble file to a list of bubbles.
     each bubble is a 4 element python list, each element is a biopython seq"""
    cortex_bubbles = []
    with open(cortex_fname) as h:
        for seq in Bio.SeqIO.parse(h, "fasta"):
            mobj = flankid_re.search(seq.id)
            if mobj:
                bubbleno = int(mobj.group(1))
                if bubbleno > len(cortex_bubbles): cortex_bubbles.append([])
                # parse the 
                if mobj.group(2) == "5p":
                    cortex_bubbles[-1].append(seq.seq)
                if mobj.group(2) == "3p":
                    cortex_bubbles[-1].append(seq.seq)
            mobj = branchid_re.search(seq.id)
            if mobj:
                cortex_bubbles[-1].append(seq.seq)
    return cortex_bubbles
            
    

sflankstr_re = re.compile(r"Start flank: (\S+) c: (\d+):(\d+)")
branchstr_re = re.compile(r"Branch: (\S+)")
eflankstr_re = re.compile(r"End flank: (\S+)")
def parse_cosmo(cosmo_fname):
    """Parse a cosmo bubble file to a list of bubbles.
     each bubble is a 4 element python list, each element is a biopython seq"""
    cosmo_bubbles = []
    with open(cosmo_fname) as h:
        for line in h:
            # check for start flank
            mobj = sflankstr_re.search(line)
            if mobj:
                cosmo_bubbles.append([Bio.Seq.Seq(mobj.group(1))])
            # check for bubble arm flank
            mobj = branchstr_re.search(line)
            if mobj:
                cosmo_bubbles[-1].append(Bio.Seq.Seq(mobj.group(1)))
            # check for end flank
            mobj = eflankstr_re.search(line)
            if mobj:
                cosmo_bubbles[-1].append(Bio.Seq.Seq(mobj.group(1)))
    return cosmo_bubbles



def diff(ref_bubbles, test_bubbles, kmer_size):
    "Check that every ref bubble occurs in the test bubble set"
    test_starts = {seq[0][-kmer_size:] for seq in test_bubbles}  

    missing_bubbles = []
    found_bubbles = []
    for bubble in ref_bubbles:
        qseq = bubble[0]
        eqseq = bubble[3]
        if str(qseq[-kmer_size:]) not in test_starts and str(eqseq[:kmer_size].reverse_complement()) not in test_starts:
            missing_bubbles.append(bubble)
        else:
            found_bubbles.append(bubble)
    return found_bubbles, missing_bubbles


# cortex_bubbles = parse_cortex("ecoli/bubbles_cortex.txt")

# cosmo_bubbles = parse_cosmo("ecoli/cosmo-color.stdout.martin")
# cosmo_bubbles2 = parse_cosmo("ecoli/cosmo-color-kmc4.stdout")
# cosmo_bubbles3 = parse_cosmo("ecoli/refactored.stdout")
cosmo_bubbles1, cosmo_bubbles2 = sys.argv[1:3]
flank_cmp_size = int(sys.argv[3])

found_bubbles, missing_bubbles = diff(parse_cosmo(cosmo_bubbles1), parse_cosmo(cosmo_bubbles2), flank_cmp_size)
print("missing bubbles:", len(missing_bubbles))
print("found bubbles:", len(found_bubbles))

