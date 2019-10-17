import subprocess
import os
import ete3
import pandas as pd
import re
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
# Load previously ran data
#
if os.path.isfile('transfers.tab'):
    transfer_df = pd.read_csv('transfers.tab', sep='\t', index_col=0)
else:
    #
    # or just do it again...
    transfers = []
    for count, group in enumerate(os.listdir('ranger')):
        if not os.path.getsize('ranger/%s/%s.output1' % (group, group)):
            #
            # if there is no RANGER's reconciliation "output1" there is no other reconciliation data, just ignore it
            continue

        #
        # otherwise, retrieve information from the reconciliation aggregate file
        transfers.append(ranger_parser.parse_aggregated(group))

    #
    # append succefully retrieved transfer data to a "master" df
    transfer_df = pd.concat(
        [tmp_transfer_df for tmp_transfer_df, tmp_gene_tree in transfers    # filter out failed parsings
                         if not tmp_transfer_df is None],                   #
        ignore_index=True,
        axis=0,
        sort=False)
#
#added this step staight to the base function, makes no sense to pass the absolute values instead of the ratio
#   it may yield a "division by zero" error, let's see...
#
#    transfer_df.loc[:, ['ranger_confidence',
#                        'ranger_confidence_donor',
#                        'ranger_confidence_recipient']] = transfer_df.loc[:, ['ranger_confidence',
#                                                                              'ranger_confidence_donor',
#                                                                              'ranger_confidence_recipient']] / 50
    #
    #save the run for future re-analysis
    transfer_df.to_csv('transfers.tab', sep='\t')

#
#filter obtained transfers based on some threshold, we are being fairly conservative here since we have a lot of data
#   there is no need to open space for spurious transfers, we can nitpick a little
#
transfer_df = transfer_df[(transfer_df.bipartition_support        >=95)  &
                          (transfer_df.ranger_confidence          >=0.8) &
                          (transfer_df.ranger_confidence_donor    >=0.8) &
                          (transfer_df.ranger_confidence_recipient>=0.8)]

########################################################################################################################
#
#preparing for the MaxTic run!
#
########################################################################################################################

#
#use a random reconciliation file as example of species tree internal branch names, as all reconciliations used the same
#   tree as model the should all be named the same
ranger_parser.name_species_tree_nodes(
    reconciliation_file='ranger/4468_BAB72593.1-GCA_000009705.1/4468_BAB72593.1-GCA_000009705.1.output1'
)

#
#write out the list of donor/recipients to pass to MaxTic
out = open('maxtic.input', 'w')
for index, row in transfer_df[['donor', 'recipient']].iterrows():
    out.write('%s\n' % '\t'.join(row.tolist()))
out.close()

ranger_parser.species_tree.write(outfile='species_tree_named_nodes',
                               format=1,
                               dist_formatter='%.10f')

#
#run MaxTic
subprocess.call(['python2.7', #needs py27... let's see for how long it will keep working...
                 '/work/ale/maxtic/MaxTiC.py',
                 'species_tree_named_nodes',
                 'maxtic.input',
                 'ls=180' #this is the value used in their paper and suggested in the documentation
                 ])

donor_recipient_pairs = []
for line in open('maxtic.input_MT_output_partial_order').readlines():
    #biggest set of donor/recipient pairs compatible to each other according to MaxTic
    donor_recipient_pairs.append('-'.join(line.split()[:2]))

#
#create a temp "donor-recipient" collumn in the df so we can use it to filter imcompatible pairs
transfer_df['donor-recipient'] = transfer_df['donor']+'-'+transfer_df['recipient']
#
#... like this!
transfer_df = transfer_df[transfer_df['donor-recipient'].isin(donor_recipient_pairs)]
#I said it was temporary, let's delete it now
transfer_df.drop(labels='donor-recipient', axis=1, inplace=True)

########################################################################################################################
#
#now that we made all the basic filtering on our transfers, let's extend the kind of data we have on them to do some
#   more fancy analysis on them
#TIP: if you are short on RAM you can just overwrite the <transfer_df>, I kept it so I can compare them in the future
extended_df = ranger_parser.assess_transfer_distance(transfer_df)
extended_df = ranger_parser.assess_dtl_cost(extended_df)
extended_df = ranger_parser.map_taxonomic_level(extended_df, taxa_table='../genomes.tab')
