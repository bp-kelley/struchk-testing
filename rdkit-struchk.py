from rdkit import RDConfig
import os, sys
from rdkit import DataStructs, Chem
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import rdStructChecker



struchk_conf_path = os.path.join(RDConfig.RDDataDir, 'struchk', '')

# don't load transform atoms for now
struchk_log_path = ''
STRUCHK_INIT = '''-tm
-tm
-or
-ca %(struchk_conf_path)scheckfgs.chk
-cc
-cl 3
-cn 999
-cs
-l %(struchk_log_path)sstruchk.log'''%locals()

opts = rdStructChecker.StructCheckerOptions()
opts.LoadGoodAugmentedAtoms("%(struchk_conf_path)scheckfgs.chk"%locals())

opts.RemoveMinorFragments = True
# -or => opts.?
opts.CheckCollisions = True
opts.CollisionLimitPercent = 3
opts.CheckStereo = True
opts.MaxMolSize = 999
opts.Verbose = False

checker = rdStructChecker.StructChecker(opts)

r = pyAvalonTools.InitializeCheckMol(STRUCHK_INIT)

def label(err):
    res = []
    for name, flag in pyAvalonTools.StruChkFlag.__dict__.items():
        try:
            if flag & err:
                res.append(name)
        except:
            pass
    return res

nitro = '''nitro.mol
  ChemDraw08311606582D

  4  3  0  0  0  0  0  0  0  0999 V2000
   -0.7145   -0.6188    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
    0.0000    0.6188    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.2062    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    0.7145   -0.6188    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  3  1  0      
  2  3  2  0      
  3  4  1  0      
M  CHG  2   1  -1   3   1
M  END
'''

i = 0
def check(ctab,f=None):
    (err, fixed_mol) = pyAvalonTools.CheckMoleculeString(ctab, False)
    mol = Chem.MolFromMolBlock(ctab, sanitize=False)
    ops = (Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_SETAROMATICITY^
           Chem.SanitizeFlags.SANITIZE_CLEANUPCHIRALITY)
    try:
        #Chem.SanitizeMol(mol, sanitizeOps=ops)
        if not mol:
            return False
    except:
        return False
    mol.UpdatePropertyCache(False)
    
    #print Chem.MolToSmiles(mol, isomericSmiles=False)
    err2 = checker.CheckMolStructure(mol)

    labels = [l.lower() for l in checker.StructureFlagsToString(err2).split(",") if l]


    if sorted(labels) == sorted(label(err)):
        print >> sys.stderr, "...ok"
        return True


    print >> sys.stderr, "...Failed" , "expected:", sorted(label(err)), "got:", sorted(labels)

    expected = repr(sorted(label(err)))
    got = repr(labels)
    ctab += "> <EXPECTED>\n%s\n\n> <GOT>\n%s\n\n$$$$\n"%(expected, got)

    oname = f
    if err:
        err = "-".join(sorted(label(err)))
        print >> sys.stderr, err
        path = os.path.join("failures", err)
        if not os.path.exists(path):
            os.mkdir(path)
        fn = os.path.join(path, oname)
        open(fn, 'w').write(ctab)
    else:
        if not os.path.exists("ok"):
            os.mkdir("ok")

        fn = os.path.join("ok", oname)
        open(fn, 'w').write(ctab)

def molcheck(fname):
    print "Examining file", fname
    text = open(fname).read()
    f = os.path.split(fname)[-1]
    mols = text.split("$$$$\n")
    print "number of molecules", len(mols)
    del text
    for i,ctab in enumerate(mols):
        check(ctab,f.replace(".", "-%06d."%i))

try:
    #check(open("bad-nitro.mol").read(), "bad-nitro-err.mol")
    #check(open("atom_check_fail.sdf").read(), "atom_check_fail-err.mol")
    #check(open("ominus.mol").read(), "ominus-err.mol")
    #check(nitro)
    #check(open("fragment.mol").read(), "fragment-err.mol")
    #check(open("fragment2.mol").read(), "fragment2-err.mol")
    #check(open("badmolblock.mol" ).read(), "badmolblock-err.mol")
    #check(open("bad_aromatic.sdf").read(), "badaromatic-err.mol")
    #check(open("atom_clash.sdf").read(), "atom_clash-err.mol")
    #molcheck("dubious.mol")
    #molcheck("atom-clash3.mol")
    #molcheck("parity.mol")
    #molcheck("dubious.sdf")
    if len(sys.argv) > 1:
        files = sys.argv[1:]
        for f in files:
            molcheck(f)
        
    else:
        for root, dirs, files in os.walk("substance"):
            for f in files:
                if ".sdf" in f:
                    print "Examining file", f
                    text = open(os.path.join(root, f)).read()
                    mols = text.split("$$$$\n")
                    print "number of molecules", len(mols)
                    del text
                    for i,ctab in enumerate(mols):
                        check(ctab,f.replace(".", "-%06d."%i))
finally:
    pyAvalonTools.CloseCheckMolFiles()
                    
                
