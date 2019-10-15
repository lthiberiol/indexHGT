import subprocess
import itertools
from copy import deepcopy
from base_functions import base_filter

#
# checks if previous run's have built a base structure needed
if not os.path.isdir('reconciliation_aggregates'):
    #
    # is there a folder with all reconciliation aggregates?
    os.mkdir('reconciliation_aggregates')

for folder in os.listdir('ranger/'):
    #
    # if there is an aggregate file for a group, keep going
    if os.path.isfile('reconciliation_aggregates/%s' % folder):
        continue

    #
    # if not copy it from its group general folder
    os.system('cp ranger/%s/aggregated reconciliation_aggregates/%s' % (folder, folder))

#
# check if there's a folder containing gene trees with already renamed tips (fixed consistencies)
if not os.path.isdir('renamed_gene_tree'):
    #
    # if not, do it...
    os.mkdir('renamed_gene_tree')

#
# check existence of such newick files, and if necessary generate them
for treefile in os.listdir('gene_trees/'):
    if not treefile.endswith('.treefile'):
        continue
    if os.path.isfile('renamed_gene_tree/%s' % treefile):
        continue

    tmp_tree = ete3.Tree('gene_trees/%s' % treefile)
    for leaf in tmp_tree.get_leaves():
        if leaf.name.count('_') == 1:
            gene, genome = leaf.name.split('_')
        elif leaf.name.count('_') > 1 and re.search('GC[AF]_', leaf.name):
            gene, genome = re.search('^([^.]+).+?(GC[AF]_\d+)', leaf.name, re.M).groups()
        elif leaf.name.count('_') == 2 and re.search('_PRJ', leaf.name):
            gene, genome = re.search('^(.+)_(PRJ.+)$', leaf.name, re.M).groups()
        else:
            print(leaf.name)
        gene = gene.split('.')[0]
        leaf.name = '%s_%s' % (genome.replace('_', ''), gene.replace('_', ''))

    tmp_tree.write(outfile='renamed_gene_tree/%s' % treefile, format=0, dist_formatter='%.10f')

#
# initiate base filtering functions
#
ranger_parser = base_filter(reference_tree              = 'species_tree-renamed',
                            gene_tree_folder            =  'renamed_gene_tree/',
                            aggregate_folder            = 'reconciliation_aggregates/',
                            reconciliation_folder       = 'ranger/',
                            overall_tree_support_thresh = 20,
                            ranger_executable           = '/work/ranger/CorePrograms/Ranger-DTL.mac',
                            leaves_allowed              = False)


#
#
#
if os.path.isfile('transfers.tab'):
    #
    transfer_df = pd.read_csv('transfers.tab', sep='\t', index_col=0)
else:
    transfers = []
    for count, group in enumerate(os.listdir('ranger')):
        if not os.path.getsize('ranger/%s/%s.output1' % (group, group)):
            continue
        transfers.append(ranger_parser.parse_aggregated(group))

    transfer_df = pd.concat([n[0] for n in transfers
                             if type(n) is not dict], ignore_index=True, axis=0, sort=False)
    transfer_df.loc[:, ['ranger_confidence',
                        'ranger_confidence_donor',
                        'ranger_confidence_recipient']] = transfer_df.loc[:, ['ranger_confidence',
                                                                              'ranger_confidence_donor',
                                                                              'ranger_confidence_recipient']] / 50

    transfer_df.to_csv('transfers.tab', sep='\t')
