import re
import sys

f = open(sys.argv[1])
fw = open(sys.argv[2],'w')
lines=f.readlines()
NewickTree=[lines[0]]

#NewickTree = ['(BMNH833953:0.16529463651919140688,(((BMNH833883:0.22945757727367316336,(BMNH724182a:0.18028180766761139897,(BMNH724182b:0.21469677818346077913,BMNH724082:0.54350916483644962085):0.00654573856803835914):0.04530853441176059537):0.02416511342888815264,(((BMNH794142:0.21236619242575086042,(BMNH743008:0.13421900772403019819,BMNH724591:0.14957653992840658219):0.02592135486124686958):0.02477670174791116522,BMNH703458a:0.22983459269245612444):0.00000328449424529074,BMNH703458b:0.29776257618061197086):0.09881729077887969892):0.02257522897558370684,BMNH833928:0.21599133163597591945):0.02365043128986757739,BMNH724053:0.16069861523756587274):0.0;']

pattern = re.compile(r"\b[0-9]+(?:\.[0-9]+)?\b")
print(NewickTree)
for tree in NewickTree:
    branch_lengths = pattern.findall(tree)
    branch_lengths.append(":")
    
    for i in range(len(branch_lengths)):
        for j in range(i+1,len(branch_lengths)):
            if branch_lengths[i] in branch_lengths[j]:
                tmp = branch_lengths[i]
                branch_lengths[i] = branch_lengths[j]
                branch_lengths[j] = tmp
    print(branch_lengths)
    new_tree=tree
    for item in branch_lengths:
        new_tree = new_tree.replace(item,"")
    fw.write(new_tree+"\n")
    print(new_tree)
fw.close()
f.close()
