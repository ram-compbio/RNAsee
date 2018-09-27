import RNA

def get_ss(rna):
    RNA.cvar.uniq_ML = 1
    subopt_data = { 'counter' : 1, 'sequence' : rna }
    # create a new model details structure
    md = RNA.md()
    # create a fold compound
    fc = RNA.fold_compound(rna, md)
    # predict Minmum Free Energy and corresponding secondary structure
    (ss, mfe) = fc.pf()
    #(ss, mfe) = fc.mfe()
    ss = fc.pbacktrack()
    # predict Minmum Free Energy and corresponding secondary structure - Partition Function
#    (ss, mfe) = fc.pf()
    # print output
    print "{}\n{} [ {} ]".format(rna, ss, mfe)
    return ss, mfe

def plot_ss(rna, ss, out_name):
    # create a new model details structure
    md = RNA.md()
    # Plot the structure to PS
    RNA.file_PS_rnaplot(rna, ss, out_name, md)
